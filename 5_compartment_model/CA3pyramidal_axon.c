/* CA3 pyramidal axon, with fixed somatic voltage.
   Bulk of model is given in Traub (1994). Compute 
   for cell n_first up to (but not including) n_last. All constants are taken 
   from Traub 1994, and then units are converted to nA, nF, uS (1/Mohm), mV, ms
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_math.h"
#include "cell.h"
#include "CA3pyramidal_axon.h"

static double C[2],g_L[2],g_Na[2],g_KDR[2];
static double V_Na = 115;
static double V_K = -25;
static double g_S;
static struct comp_conn c_conn[5][2];

static struct CA3pyramidal_axon rhs(struct CA3pyramidal_axon* axon,double t);
static double a_m(double V);
static double b_m(double V);
static double a_h(double V);
static double b_h(double V);
static double a_n(double V);
static double b_n(double V);
static int level(int comp);
static void compartment_connections(struct comp_conn[][2]);

void CA3pyraxon_setup_constants(){
  /* Setup constants for pyramical cell*/

  // allocate variables and initialize some constants
  // set up numbers by level
  double A[2] = {942,236};
  double Na = 500;
  double KDR = 250;
  double Rm_A = 1000; //membrane resistance for axon and IS
  extern double C[2],g_L[2],g_Na[2],g_KDR[2];
  extern struct comp_conn c_conn[5][2];
  
  int k;  

  // initialize rest of constants
  for(k=0;k<2;k++){
    C[k] = .75e-5*A[k];
    g_Na[k] = 1e-5*Na*A[k];
    g_KDR[k] = 1e-5*KDR*A[k];
    g_L[k] = 1e-2*A[k]/Rm_A; 
  }
  
  // sets conductance from soma to axon as well
  compartment_connections(c_conn);
}

void CA3pyraxon_init(struct CA3pyramidal_axon* axon, int id, double V_S){
  // Initialize the pyramical cell pointed to by pyr 
  int k;

  axon->id = id;

  // set initial data for neuron
  // resting V,m,h,n for single neuron with no input, assume V=0
  for(k=0;k<5;k++){
    axon->V[k] = 0;
    axon->m[k] = a_m(axon->V[k])/(a_m(axon->V[k])+b_m(axon->V[k]));
    axon->h[k] = a_h(axon->V[k])/(a_h(axon->V[k])+b_h(axon->V[k]));
    axon->n[k] = a_n(axon->V[k])/(a_n(axon->V[k])+b_n(axon->V[k]));
  }

  //printf("V_S = %g\n",V_S);
  axon->V_S = V_S;

  axon->currents = NULL;
  axon->n_currents = 0;
}

void CA3pyraxon_set_currents(struct CA3pyramidal_axon* axon, 
			  struct current *currents, 
			  int n_currents){
  axon->currents = currents;
  axon->n_currents = n_currents;
}

void CA3pyraxon_step(struct CA3pyramidal_axon* axon0, 
		     struct CA3pyramidal_axon* axon1, 
		     double dt, double t){
  // take axon0 and integrate using midpoint method for one step, 
  // putting the answer in axon1
  struct CA3pyramidal_axon axon_h; // pyramidal half step state
  struct CA3pyramidal_axon der; // derivative for all variables
  int k;

  // then apply midpoint method
  //printf("Calling rhs\n");
  der = rhs(axon0,t);
  //printf("computing half step\n");
  for(k=0;k<5;k++){
    axon_h.m[k] = .5*dt*der.m[k] + axon0->m[k];
    axon_h.h[k] = .5*dt*der.h[k] + axon0->h[k];
    axon_h.n[k] = .5*dt*der.n[k] + axon0->n[k];
    axon_h.V[k] = .5*dt*der.V[k] + axon0->V[k];
  }

  //printf("calling rhs again\n");
  axon_h.id = axon0->id;
  axon_h.n_currents = axon0->n_currents;
  axon_h.currents = axon0->currents;
  axon_h.V_S = axon0->V_S;
  der = rhs(&axon_h,t+dt/2);
  //printf("computing axon1\n");
  for(k=0;k<5;k++){
    axon1->m[k] = dt*der.m[k] + axon0->m[k];
    axon1->h[k] = dt*der.h[k] + axon0->h[k];
    axon1->n[k] = dt*der.n[k] + axon0->n[k];
    axon1->V[k] = dt*der.V[k] + axon0->V[k];
  }
    
  return;
}
 
struct CA3pyramidal_axon rhs(struct CA3pyramidal_axon* axon,double t){
   /* Evaluate the right-hand-side of the differential equation
      that governs the evolution of a CA3pyramidal_axon */
  extern double C[2],g_L[2],g_Na[2],g_KDR[2];
  extern struct comp_conn c_conn[5][2];
  extern double g_S;
  struct CA3pyramidal_axon der; // derivative/right-hand-side
  int j,k;
  double I;

  //printf("calculating for IS and axon\n");
  //printf("axon->n_currents = %d\n",axon->n_currents);
  // then calculate for IS and axon
  for(k=0;k<5;k++){
    der.m[k] = a_m(axon->V[k])*(1-axon->m[k])
    - b_m(axon->V[k])*axon->m[k];
    der.h[k] = a_h(axon->V[k])*(1-axon->h[k])
      - b_h(axon->V[k])*axon->h[k];
    der.n[k] = a_n(axon->V[k])*(1-axon->n[k])
      - b_n(axon->V[k])*axon->n[k];
  }
  der.V[0] = g_L[0]*(0-axon->V[0])
      + g_Na[0]*pow(axon->m[0],3)*axon->h[0]*(V_Na-axon->V[0])
      + g_KDR[0]*pow(axon->n[0],4)*(V_K-axon->V[0]);
  for(k=1;k<5;k++){ 
   der.V[k] = g_L[1]*(0-axon->V[k])
      + g_Na[1]*pow(axon->m[k],3)*axon->h[k]*(V_Na-axon->V[k])
      + g_KDR[1]*pow(axon->n[k],4)*(V_K-axon->V[k]);
  }

  for(k=0;k<axon->n_currents;k++){
    // only call cur_func once for Poisson process
    I = (*axon->currents[k].cur_func)(axon->id,t);
    der.V[axon->currents[k].comp] = der.V[axon->currents[k].comp] + I;
      //(*axon->currents[k].cur_func)(axon->id,t);
  }

  der.V[0] = der.V[0] + g_S*(axon->V_S-axon->V[0]);
  for(k=0;k<5;k++){
    j = 0;
    while(j<2 && c_conn[k][j].comp != -1){
      der.V[k] = der.V[k] + c_conn[k][j].gamma*
	(axon->V[c_conn[k][j].comp]-axon->V[k]);
      j = j+1;
    }
  }
  der.V[0] = der.V[0]/C[0];
  for(k=1;k<5;k++){
    der.V[k] = der.V[k]/C[1];
  }
  
  return der;
}

// functions for the Initial Segment and Axonal compartments
double a_m(double V){
  double ans = .8*(17.2-V)/(exp((17.2-V)/4) - 1);
  return ans;
}
double b_m(double V){
  double ans = .7*(V-42.2)/(exp((V-42.2)/5) - 1);
  return ans;
}
double a_h(double V){
  double ans = .32*exp((42-V)/18);
  return ans;
}
double b_h(double V){
  double ans = 10/(exp((42-V)/5) + 1);
  return ans;
}
double a_n(double V){
  double ans = .03*(17.2-V)/(exp((17.2-V)/5) - 1);
  return ans;
}
double b_n(double V){
  double ans = .45*exp((12-V)/40);
  return ans;
}

// Rows of matrix correspond to compartment number, and hold
// structs describing connections to other compartments.
void compartment_connections(struct comp_conn c_conn[][2]){
  // r=radius in micro-meters
  double r[] = {15,2,.5};
  // L=length in micro-meters
  double L[] = {25.5,75,75};
  // R= internal resistance in Ohm-cm
  double R[3];
  double rho[3];
  double g[2];
  extern double g_S;
  int i,j;

  R[0] = 200;
  R[1] = 100;
  R[2] = 100;
  
  for(i=0;i<3;i++){
    rho[i] = R[i]*L[i]*1e-2/(PI*pow(r[i],2));
  }

  g_S = 2/(rho[0]+rho[1]); //conductance from soma to IS
  //printf("g_S = %g\n",g_S);
  g[0] = 2/(rho[1]+rho[2]); //conductance from IS to rest of axon
  g[1] = 1/rho[2]; //conductance between axonal compartments
  
  // first initialize array
  for(i=0;i<5;i++){
    for(j=0;j<2;j++){
      c_conn[i][j] = make_comp_conn(-1,0);
    }
  }
  // then make compartment connections
  c_conn[0][0] = make_comp_conn(1,g[0]);
  c_conn[1][0] = make_comp_conn(0,g[0]);
  c_conn[1][1] = make_comp_conn(2,g[1]);
  c_conn[2][0] = make_comp_conn(1,g[1]);
  c_conn[2][1] = make_comp_conn(3,g[1]);
  c_conn[3][0] = make_comp_conn(2,g[1]);
  c_conn[3][1] = make_comp_conn(4,g[1]);
  c_conn[4][0] = make_comp_conn(3,g[1]);
  return;
}
