/* give the voltage of a network of neurons in mV
using 69 compartments (soma, initial segment and axon) based on Traub's model
(1999) for the given run_time using time step dt. The network has "cells"
number of cells in an n_row by
n_col grid with connections conn. All connections
are modeled as gap junctions. Bulk of model is given
in Traub (1994). Compute for cell n_first up to (but not including) n_last.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979

struct comp_conn {
  int comp;
  double gamma;
};

double a_m_S(double V);
double b_m_S(double V);
double a_h_S(double V);
double b_h_S(double V);
double a_s_S(double V);
double b_s_S(double V);
double a_n_S(double V);
double b_n_S(double V);
double a_q_S(double X);
double b_q_S();
double a_a_S(double V);
double b_a_S(double V);
double a_b_S(double V);
double b_b_S(double V);
double a_c_S(double V);
double b_c_S(double V);
double a_m(double V);
double b_m(double V);
double a_h(double V);
double b_h(double V);
double a_n(double V);
double b_n(double V);
int level(int comp);
void compartment_connections(struct comp_conn[][6]);
struct comp_conn make_comp_conn(int comp,double gamma);
double min(double a,double b);
int my_round(double x);

void traub69(int cells,int conn[][4],double run_time,int start_rec, 
	     double dt,int* cell_assgn,int rank,int nRanks){
  /* All constants are taken from Traub 1999, and then units are 
     converted to nA, nF, uS (1/Mohm), mV, ms
  */
  
  // allocate variables and initialize some constants
  // set up numbers by level
  double A[15] = {880,792,471,2403,3142,1979,1571,1188,942,603,603,
		  1223,440,942,236};
  double Na[15] = {0,0,1,100,3,3,0,0,0,0,0,0,0,500,500};
  double Ca[13] = {1,1,1,1,1,1,2,3,3,1,1,2,3};
  double KDR[15] = {0,0,15,135,20,20,0,0,0,0,0,0,0,250,250};
  double KAHP = .8; 
  double KC[13] = {4,4,8,20,8,8,4,12,12,4,4,4,12};
  double KA[13] = {0,.5,.5,.5,.5,.5,0,0,0,0,0,0,0};
  double phi_init[13] = {148,164,123,24,18,29,37,13,16,25,25,47,34};
  double C[69],g_L[69],g_Na[69],g_KDR[69];
  double g_Ca[64],g_KAHP[64],g_KC[64],g_KA[64],phi[64];
  double Rm_SD = 50000; // membrane resistance for soma and dendrites
  double Rm_A = 1000; //membrane resistance for axon and IS
  double beta_S = .001; // inverse time constant for chi (X) for soma
  double beta_D = .05; // inverse time constant for chi (X) for dendrites
  double g_gj = 3.66e-3;
  double V_Na = 115;
  double V_Ca = 140;
  double V_K_S = -15;
  double V_K = -25;
  double lambda = .002; // stimulus/ms/cell %1
  double t_stim = .3125; // amount of time to stimulate %.25
  int rec_step,stim_step;
    
  int n_first = cell_assgn[rank];
  int n_last = cell_assgn[rank+1];
  int n_total = n_last - n_first;

  // Only keep track of everything for this node.
  double V[3][n_total][69];
  double m[3][n_total][69];
  double h[3][n_total][69];  
  double n[3][n_total][69];
  // variables for soma/dendrites
  double s[3][n_total][64];
  double a[3][n_total][64];
  double b[3][n_total][64];
  double q[3][n_total][64];
  double c[3][n_total][64],X[3][n_total][64];
   // Keep track of penultimate axonal compartment for everyone
  double Vpac[cells];
  double t_count[n_total];
  double I_soma[n_total];
  struct comp_conn c_conn[69][6];
  
  MPI_Status status;
  char file[100];
  int x,k,i,j,lev;
  
  FILE *outp; /* pointer to output file */
  FILE *outstim;

  /* Setup output file */
  sprintf(file,"/scratch/emunro01/t69_%d_%d-%d.out",rank,my_round(start_rec),
	  my_round(start_rec)+499);
  outp = fopen(file,"w");
  sprintf(file,"/scratch/emunro01/t69stim_%d_%d-%d.out",rank,
	  my_round(start_rec),my_round(start_rec)+499);
  outstim = fopen(file,"w");
  printf("%d: Setup output file.\n",rank);
  fflush(NULL);

  /* Figure out when to record data: rec_step */
  rec_step = my_round(.1*2/dt);
  stim_step = my_round(100*2/dt);

  // initialize rest of constants
  for(k=0;k<69;k++){
    lev = level(k);
    C[k] = .75e-5*A[lev];
    g_Na[k] = 1e-5*Na[lev]*A[lev];
    g_KDR[k] = 1e-5*KDR[lev]*A[lev];
  }
  
  for(k=0;k<64;k++){
    lev = level(k);
    // leak conductance is 1/Rm by Traub's code, lines 1667 and 2205
    g_L[k] = 1e-2*A[lev]/Rm_SD; 
    g_Ca[k] = 1e-5*Ca[lev]*A[lev]/2; // Ca conductance is halved for T1999
    g_KAHP[k] = 1e-5*KAHP*A[lev];
    g_KC[k] = 1e-5*KC[lev]*A[lev];
    g_KA[k] = 1e-5*KA[lev]*A[lev];
    phi[k] = phi_init[lev];
  }
	
  for(k=64;k<69;k++){
    lev = level(k);
    g_L[k] = 1e-2*A[lev]/Rm_A; 
  }
	
  // set initial data for neuron
  // resting V,m,h,n for single neuron with no input, assume V=0
  V[0][0][0] = 0; 
  m[0][0][0] = a_m_S(V[0][0][0])/(a_m_S(V[0][0][0])+b_m_S(V[0][0][0]));
  h[0][0][0] = a_h_S(V[0][0][0])/(a_h_S(V[0][0][0])+b_h_S(V[0][0][0]));
  n[0][0][0] = a_n_S(V[0][0][0])/(a_n_S(V[0][0][0])+b_n_S(V[0][0][0]));
  s[0][0][0] = a_s_S(V[0][0][0])/(a_s_S(V[0][0][0])+b_s_S(V[0][0][0]));
  X[0][0][0] = phi[0]*g_Ca[0]*pow(s[0][0][0],2)*(V_Ca-V[0][0][0])/beta_D;
  q[0][0][0] = a_q_S(X[0][0][0])/(a_q_S(X[0][0][0])+b_q_S());
  a[0][0][0] = a_a_S(V[0][0][0])/(a_a_S(V[0][0][0])+b_a_S(V[0][0][0]));
  b[0][0][0] = a_b_S(V[0][0][0])/(a_b_S(V[0][0][0])+b_b_S(V[0][0][0]));
  c[0][0][0] = a_c_S(V[0][0][0])/(a_c_S(V[0][0][0])+b_c_S(V[0][0][0]));
  V[0][0][64] = V[0][0][0];
  m[0][0][64] = a_m(V[0][0][64])/(a_m(V[0][0][64])+b_m(V[0][0][64]));
  h[0][0][64] = a_h(V[0][0][64])/(a_h(V[0][0][64])+b_h(V[0][0][64]));
  n[0][0][64] = a_n(V[0][0][64])/(a_n(V[0][0][64])+b_n(V[0][0][64]));
  for(x=0;x<n_total;x++){
    for(i=0;i<=1;i++){
      for(k=0;k<64;k++){
        V[i][x][k] = V[0][0][0];
        m[i][x][k] = m[0][0][0];
        h[i][x][k] = h[0][0][0];
        n[i][x][k] = n[0][0][0];
        s[i][x][k] = s[0][0][0];
        X[i][x][k] = phi[k]*g_Ca[k]*pow(s[0][0][0],2)*(V_Ca-V[0][0][0])/beta_D;
        q[i][x][k] = q[0][0][0];
        a[i][x][k] = a[0][0][0];
        b[i][x][k] = b[0][0][0];
        c[i][x][k] = c[0][0][0];
      }
      X[i][x][28] = phi[28]*g_Ca[28]*pow(s[0][0][0],2)*(V_Ca-V[0][0][0])
	/beta_S;
      for(k=64;k<69;k++){
        V[i][x][k] = V[0][0][64];
        m[i][x][k] = m[0][0][64];
        h[i][x][k] = h[0][0][64];
        n[i][x][k] = n[0][0][64];
      }
    }
  }

  for(x=0;x<cells;x++){
    Vpac[x] = V[0][0][0];
  }
  for(x=0;x<n_total;x++){
    t_count[x] = 0;
    srand(x+n_first);
    I_soma[x] = .05 + .01*rand()/RAND_MAX;
  }
  compartment_connections(c_conn);

  printf("%d: proceeding with midpoint method.\n",rank);
  fflush(NULL);
  // proceed with midpoint method
  for(i=1;i<=ceil(2*run_time/dt);i++){
    if(fmod(i,rec_step) == 0 && i>=my_round(start_rec*2/dt)){
      if(fmod(i,my_round(500*2/dt)) == 0 && i!=my_round(start_rec*2/dt) 
         &&i!=ceil(2*run_time/dt)){
        // Setup output file 
        fclose(outp);
        fclose(outstim);
        sprintf(file,"/scratch/emunro01/t69_%d_%d-%d.out",rank,my_round(i*dt/2),
                my_round(i*dt/2)+499);
        outp = fopen(file,"w");
        sprintf(file,"/scratch/emunro01/t69stim_%d_%d-%d.out",rank,my_round(i*dt/2),
                my_round(i*dt/2)+499);
        outstim = fopen(file,"w");
      }
      fprintf(outp,"%g",dt*i/2);
    }
    for(x=0;x<n_total;x++){
      srand(cells*i + x+n_first);
      // calculate for soma/dendrites first
      // calculate right side of diffeq
      for(k=0;k<64;k++){
	m[2][x][k] = a_m_S(V[1][x][k])*(1-m[1][x][k])
	  - b_m_S(V[1][x][k])*m[1][x][k];
	h[2][x][k] = a_h_S(V[1][x][k])*(1-h[1][x][k])
	  - b_h_S(V[1][x][k])*h[1][x][k];
	s[2][x][k] = a_s_S(V[1][x][k])*(1-s[1][x][k])
	  - b_s_S(V[1][x][k])*s[1][x][k];
	n[2][x][k] = a_n_S(V[1][x][k])*(1-n[1][x][k])
	  - b_n_S(V[1][x][k])*n[1][x][k];
	if(k==28){
	  X[2][x][k] = phi[k]*g_Ca[k]*pow(s[1][x][k],2)*(V_Ca-V[1][x][k])
	    - beta_S*X[1][x][k];
	}else{
	  X[2][x][k] = phi[k]*g_Ca[k]*pow(s[1][x][k],2)*(V_Ca-V[1][x][k])
	    - beta_D*X[1][x][k];
	}
	q[2][x][k] = a_q_S(X[1][x][k])*(1-q[1][x][k])
	  - b_q_S(X[1][x][k])*q[1][x][k];
	a[2][x][k] = a_a_S(V[1][x][k])*(1-a[1][x][k])
	  - b_a_S(V[1][x][k])*a[1][x][k];
	b[2][x][k] = a_b_S(V[1][x][k])*(1-b[1][x][k])
	  - b_b_S(V[1][x][k])*b[1][x][k];
	c[2][x][k] = a_c_S(V[1][x][k])*(1-c[1][x][k])
	  - b_c_S(V[1][x][k])*c[1][x][k];
	V[2][x][k] = g_L[k]*(0-V[1][x][k])
	  + g_Na[k]*pow(m[1][x][k],2)*h[1][x][k]*(V_Na-V[1][x][k])
	  + g_Ca[k]*pow(s[1][x][k],2)*(V_Ca-V[1][x][k])
	  + g_KDR[k]*pow(n[1][x][k],2)*(V_K_S-V[1][x][k])
	  + g_KA[k]*a[1][x][k]*b[1][x][k]*(V_K_S-V[1][x][k])
	  + g_KAHP[k]*q[1][x][k]*(V_K_S-V[1][x][k])
	  + g_KC[k]*c[1][x][k]*min(1,X[1][x][k]/250)*(V_K_S-V[1][x][k]);
	if(k==28){
	  V[2][x][k] = V[2][x][k] + I_soma[x]; 
	  // dep current according to T1999
	}
	j = 0;
	while(j<6 && c_conn[k][j].comp != -1){
	  V[2][x][k] = V[2][x][k] + c_conn[k][j].gamma*
	    (V[1][x][c_conn[k][j].comp]-V[1][x][k]);
	  j = j+1;
	}
	V[2][x][k] = V[2][x][k]/C[k];
	
	// then apply midpoint method
	if(fmod(i,2) == 0){
	  m[2][x][k] = dt*m[2][x][k] + m[0][x][k];
	  h[2][x][k] = dt*h[2][x][k] + h[0][x][k];
          s[2][x][k] = dt*s[2][x][k] + s[0][x][k];
          n[2][x][k] = dt*n[2][x][k] + n[0][x][k];
          X[2][x][k] = dt*X[2][x][k] + X[0][x][k];
          q[2][x][k] = dt*q[2][x][k] + q[0][x][k];
          a[2][x][k] = dt*a[2][x][k] + a[0][x][k];
          b[2][x][k] = dt*b[2][x][k] + b[0][x][k];
          c[2][x][k] = dt*c[2][x][k] + c[0][x][k];
          V[2][x][k] = dt*V[2][x][k] + V[0][x][k];
        }else{
          m[2][x][k] = .5*dt*m[2][x][k] + m[1][x][k];
          h[2][x][k] = .5*dt*h[2][x][k] + h[1][x][k];
          s[2][x][k] = .5*dt*s[2][x][k] + s[1][x][k];
          n[2][x][k] = .5*dt*n[2][x][k] + n[1][x][k];
          X[2][x][k] = .5*dt*X[2][x][k] + X[1][x][k];
          q[2][x][k] = .5*dt*q[2][x][k] + q[1][x][k];
          a[2][x][k] = .5*dt*a[2][x][k] + a[1][x][k];
          b[2][x][k] = .5*dt*b[2][x][k] + b[1][x][k];
          c[2][x][k] = .5*dt*c[2][x][k] + c[1][x][k];
          V[2][x][k] = .5*dt*V[2][x][k] + V[1][x][k];
        }
      }
      
      // then calculate for IS and axon
      for(k=64;k<69;k++){
        m[2][x][k] = a_m(V[1][x][k])*(1-m[1][x][k])
	  - b_m(V[1][x][k])*m[1][x][k];
	h[2][x][k] = a_h(V[1][x][k])*(1-h[1][x][k])
	  - b_h(V[1][x][k])*h[1][x][k];
	n[2][x][k] = a_n(V[1][x][k])*(1-n[1][x][k])
	  - b_n(V[1][x][k])*n[1][x][k];
	V[2][x][k] = g_L[k]*(0-V[1][x][k])
          + g_Na[k]*pow(m[1][x][k],3)*h[1][x][k]*(V_Na-V[1][x][k])
          + g_KDR[k]*pow(n[1][x][k],4)*(V_K-V[1][x][k]);
        if(k == 67){
          // figure in current from gap junctions
          j = 0;
          while(j<4 && conn[x+n_first][j] != -1){
            V[2][x][k] = V[2][x][k] + g_gj*
            (Vpac[conn[x+n_first][j]]-V[1][x][k]);
            j = j+1;
          }
        }

        if(k == 68){
          // see about spontaneous stimulation
          if((double)rand()/RAND_MAX < (dt/2)*lambda*exp(-(dt/2)*lambda)){
            // Poisson probability of stimulus during dt
            // given lambda in stim/ms/cell.
            if(i >= start_rec*2/dt){
	      fprintf(outstim,"%d\t%g\n",x+n_first,dt*i/2);
	      fflush(NULL);
	    }
            V[2][x][k] = V[2][x][k] + .2;
            t_count[x] = t_stim;
          }else if(t_count[x] > 0){
            V[2][x][k] = V[2][x][k] + .2;
            t_count[x] = t_count[x] - dt/2;
          }
        }

        j = 0;
        while(j<6 && c_conn[k][j].comp != -1){
          V[2][x][k] = V[2][x][k] + c_conn[k][j].gamma*
          (V[1][x][c_conn[k][j].comp]-V[1][x][k]);
          j = j+1;
        }
        V[2][x][k] = V[2][x][k]/C[k];

        if(fmod(i,2) == 0){
          m[2][x][k] = dt*m[2][x][k] + m[0][x][k];
          h[2][x][k] = dt*h[2][x][k] + h[0][x][k];
          n[2][x][k] = dt*n[2][x][k] + n[0][x][k];
          V[2][x][k] = dt*V[2][x][k] + V[0][x][k];
        }else{
          m[2][x][k] = .5*dt*m[2][x][k] + m[1][x][k];
          h[2][x][k] = .5*dt*h[2][x][k] + h[1][x][k];
          n[2][x][k] = .5*dt*n[2][x][k] + n[1][x][k];
          V[2][x][k] = .5*dt*V[2][x][k] + V[1][x][k];
        }
      } // for k=1;k<6
    } // for x=0;x<n_total
    
    /* Update voltages after whole network is calculated to get
    networking right. */
  
    for(x=0;x<n_total;x++){
      Vpac[x+n_first] = V[2][x][67];
      if(fmod(i,rec_step) == 0 && i>=start_rec*2/dt){
        fprintf(outp,"\t%g",V[2][x][28]);
        fprintf(outp,"\t%g",V[2][x][67]);
        fflush(outp);
      }
    }

    // Send/Receive voltages needed for gap junctions.
    for(x=0;x<rank;x++){
      MPI_Send(&Vpac[n_first],n_total,MPI_DOUBLE,x,1,MPI_COMM_WORLD);
    }
    for(x=rank+1;x<nRanks;x++){
      MPI_Send(&Vpac[n_first],n_total,MPI_DOUBLE,x,1,MPI_COMM_WORLD);
    }
    for(x=0;x<rank;x++){
      MPI_Recv(&Vpac[cell_assgn[x]],cell_assgn[x+1]-cell_assgn[x],
	       MPI_DOUBLE,x,1,MPI_COMM_WORLD,&status);
    }
    for(x=rank+1;x<nRanks;x++){
      MPI_Recv(&Vpac[cell_assgn[x]],cell_assgn[x+1]-cell_assgn[x],
	       MPI_DOUBLE,x,1,MPI_COMM_WORLD,&status);
    }

    memcpy(&V[1][0][0],&V[2][0][0],sizeof(double)*n_total*69);
    memcpy(&m[1][0][0],&m[2][0][0],sizeof(double)*n_total*69);
    memcpy(&n[1][0][0],&n[2][0][0],sizeof(double)*n_total*69);
    memcpy(&h[1][0][0],&h[2][0][0],sizeof(double)*n_total*69);
    memcpy(&s[1][0][0],&s[2][0][0],sizeof(double)*n_total*64);
    memcpy(&X[1][0][0],&X[2][0][0],sizeof(double)*n_total*64);
    memcpy(&q[1][0][0],&q[2][0][0],sizeof(double)*n_total*64);
    memcpy(&a[1][0][0],&a[2][0][0],sizeof(double)*n_total*64);
    memcpy(&b[1][0][0],&b[2][0][0],sizeof(double)*n_total*64);
    memcpy(&c[1][0][0],&c[2][0][0],sizeof(double)*n_total*64);
    if(fmod(i,2)==0){
      memcpy(&V[0][0][0],&V[2][0][0],sizeof(double)*n_total*69);
      memcpy(&m[0][0][0],&m[2][0][0],sizeof(double)*n_total*69);
      memcpy(&n[0][0][0],&n[2][0][0],sizeof(double)*n_total*69);
      memcpy(&h[0][0][0],&h[2][0][0],sizeof(double)*n_total*69);
      memcpy(&s[0][0][0],&s[2][0][0],sizeof(double)*n_total*64);
      memcpy(&X[0][0][0],&X[2][0][0],sizeof(double)*n_total*64);
      memcpy(&q[0][0][0],&q[2][0][0],sizeof(double)*n_total*64);
      memcpy(&a[0][0][0],&a[2][0][0],sizeof(double)*n_total*64);
      memcpy(&b[0][0][0],&b[2][0][0],sizeof(double)*n_total*64);
      memcpy(&c[0][0][0],&c[2][0][0],sizeof(double)*n_total*64);
    }

    /* Start new line for next time step */
    if(fmod(i,rec_step) == 0 && i>=start_rec*2/dt){
      fprintf(outp,"\n");
    }
    
  } // for i = 1:lg
  
  fclose(outp);
  fclose(outstim);
  return;
}

// functions for the Soma
double a_m_S(double V){
  double ans = .32*(13.1-V)/(exp((13.1-V)/4)-1);
  return ans;
}
double b_m_S(double V){
  double ans = .28*(V-40.1)/(exp((V-40.1)/5) - 1);
  return ans;
}
double a_h_S(double V){
  double ans = .128*exp((17-V)/18);
  return ans;
}
double b_h_S(double V){
  double ans = 4/(1+exp((40-V)/5));
  return ans;
}
double a_s_S(double V){
  double ans = 1.6/(1 + exp(-.072*(V-65)));
  return ans;
}
double b_s_S(double V){
  double ans =.02*(V-51.1)/(exp((V-51.1)/5)-1);
  return ans;
}
double a_n_S(double V){
  double ans = .016*(35.1-V)/(exp((35.1-V)/5) - 1);
  return ans;
}
double b_n_S(double V){
  double ans = .25*exp((20-V)/40);
  return ans;
}
double a_q_S(double X){
  double ans = min(.2e-4*X,.01);
  return ans;
}
double b_q_S(){
  double ans = .001;
  return ans;
}
double a_a_S(double V){
  double ans = .02*(13.1-V)/(exp((13.1-V)/10) - 1);
  return ans;
}
double b_a_S(double V){
  double ans = .0175*(V-40.1)/(exp((V-40.1)/10) - 1);
  return ans;
}
double a_b_S(double V){
  double ans = .0016*exp((-13-V)/18);
  return ans;
}
double b_b_S(double V){
  double ans = .05/(1 + exp((10.1-V)/5));
  return ans;
}
double a_c_S(double V){
  double ans;
  if(V <= 50){
    ans = exp((V-10)/11 - (V-6.5)/27)/18.975;
  }else{
    ans = 2*exp((6.5-V)/27);
  }
  return ans;
}
double b_c_S(double V){
  double ans;
  if(V <= 50){
    ans = 2*exp((6.5-V)/27) - a_c_S(V);
  }else{
    ans = 0;
  }
  return ans;
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

// Helper function to map compartment to level
int level(int comp){
  int lev;
  if(comp <= 15){
    lev = 0; //level 1: level-1
  }else if(comp>=16 && comp<=23){
    lev = 1;
  }else if(comp>=24 && comp<=27){
    lev = 2;
  }else if(comp==28){
    lev = 3; //Soma
  }else if(comp==29){
    lev = 4;
  }else if(comp==30 || comp==31){
    lev = 5;
  }else if(comp==33 || comp==34){
    lev = 6;
  }else if(comp>=36 && comp<=39){
    lev = 7;
  }else if(comp==41 || comp==42 || comp==45 || comp==46){
    lev = 8;
  }else if(comp>=48 && comp<=55){
    lev = 9;
  }else if(comp>=56 && comp<=63){
    lev = 10;
  }else if(comp==32 || comp==35){
    lev = 11; //7 oblique
  }else if(comp==40 || comp==43 || comp==44 || comp==47){
    lev = 12; //9 oblique
  }else if(comp==64){
    lev = 13; //IS
  }else if(comp>=65 && comp<=68){
    lev = 14; //axon
  }
  return lev;
}

// Rows of matrix correspond to compartment number, and hold
// structs describing connections to other compartments.
void compartment_connections(struct comp_conn c_conn[][6]){
  // r=radius in micro-meters
  double r[] = {1,1.57,2.5,15,5,3.15,2.5,1.57,1.25,.8,.8,1.39,.7,2,.5};
  // L=length in micro-meters
  double L[] = {70,40,15,25.5,50,50,50,60,60,60,60,70,50,75,75};
  // R= internal resistance in Ohn-cm
  double R[15];
  double rho[15];
  double g[27];
  int i,j;

  for(i=0;i<13;i++){
    R[i] = 200;
  }
  R[13] = 100;
  R[14] = 100;
  for(i=0;i<15;i++){
    rho[i] = R[i]*L[i]*1e-2/(PI*pow(r[i],2));
  }
  for(i=0;i<10;i++){ //uS
    g[i] = 1/rho[i]; //lev i to itself
    g[i+10] = 2/(rho[i]+rho[i+1]); //lev i to i+1
  }
  g[20] = 2/(rho[5] + rho[11]); // 6-7 oblique
  g[21] = 2/(rho[6] + rho[11]); // 7-7 oblique
  g[22] = 2/(rho[7] + rho[12]); // 8-9 oblique
  g[23] = 2/(rho[8] + rho[12]); // 9-9 oblique
  g[24] = 2/(rho[3] + rho[13]); // soma-IS
  g[25] = 2/(rho[13] + rho[14]); // IS-axon
  g[26] = 1/rho[14]; // axon-axon;

  // first initialize array
  for(i=0;i<69;i++){
    for(j=0;j<6;j++){
      c_conn[i][j] = make_comp_conn(-1,0);
    }
  }
  // then set individual compartment connections
  for(i=0;i<15;i=i+2){ //level 1
    c_conn[i][0] = make_comp_conn(i+1,g[0]);
    c_conn[i][1] = make_comp_conn(16+i/2,g[10]);
    c_conn[i+1][0] = make_comp_conn(i,g[0]);
    c_conn[i+1][1] = make_comp_conn(16+i/2,g[10]);
  }
  for(i=0;i<8;i=i+2){ //level 2
    c_conn[16+i][0] = make_comp_conn(i*2,g[10]);
    c_conn[16+i][1] = make_comp_conn(i*2+1,g[10]);
    c_conn[16+i][2] = make_comp_conn(16+i+1,g[1]);
    c_conn[16+i][3] = make_comp_conn(24+i/2,g[11]);
    c_conn[16+i+1][0] = make_comp_conn((i+1)*2,g[10]);
    c_conn[16+i+1][1] = make_comp_conn((i+1)*2+1,g[10]);
    c_conn[16+i+1][2] = make_comp_conn(16+i,g[1]);
    c_conn[16+i+1][3] = make_comp_conn(24+i/2,g[11]);
  }
  for(i=0;i<4;i++){ //level 3, assume each attached to soma separately
    c_conn[24+i][0] = make_comp_conn(16+i*2,g[11]);
    c_conn[24+i][1] = make_comp_conn(16+i*2+1,g[11]);
    c_conn[24+i][2] = make_comp_conn(28,g[12]);
  }
  for(i=0;i<4;i++){//level 4 to level 5
    c_conn[28][i] = make_comp_conn(24+i,g[12]);
  }
  c_conn[28][4] = make_comp_conn(29,g[13]);
  c_conn[28][5] = make_comp_conn(64,g[24]);
  c_conn[29][0] = make_comp_conn(28,g[13]);
  c_conn[29][1] = make_comp_conn(30,g[14]);
  c_conn[29][2] = make_comp_conn(31,g[14]);
  c_conn[30][0] = make_comp_conn(29,g[14]);
  c_conn[30][1] = make_comp_conn(32,g[20]);
  c_conn[30][2] = make_comp_conn(33,g[15]);
  c_conn[31][0] = make_comp_conn(29,g[14]);
  c_conn[31][1] = make_comp_conn(34,g[15]);
  c_conn[31][2] = make_comp_conn(35,g[20]);
  c_conn[32][0] = make_comp_conn(30,g[20]);
  c_conn[32][1] = make_comp_conn(33,g[21]);
  c_conn[33][0] = make_comp_conn(30,g[15]);
  c_conn[33][1] = make_comp_conn(32,g[21]);
  c_conn[33][2] = make_comp_conn(36,g[16]);
  c_conn[33][3] = make_comp_conn(37,g[16]);
  c_conn[34][0] = make_comp_conn(31,g[15]);
  c_conn[34][1] = make_comp_conn(35,g[21]);
  c_conn[34][2] = make_comp_conn(38,g[16]);
  c_conn[34][3] = make_comp_conn(39,g[16]);
  c_conn[35][0] = make_comp_conn(31,g[20]);
  c_conn[35][1] = make_comp_conn(34,g[21]);
  for(i=0;i<4;i=i+2){ //level 8
    c_conn[36+i][0] = make_comp_conn(33+i/2,g[16]);
    c_conn[36+i][1] = make_comp_conn(36+i+1,g[7]);
    c_conn[36+i][2] = make_comp_conn(40+i*2,g[22]);
    c_conn[36+i][3] = make_comp_conn(40+i*2+1,g[17]);
    c_conn[36+i+1][0] = make_comp_conn(33+i/2,g[16]);
    c_conn[36+i+1][1] = make_comp_conn(36+i,g[7]);
    c_conn[36+i+1][2] = make_comp_conn(40+(i+1)*2,g[17]);
    c_conn[36+i+1][3] = make_comp_conn(40+(i+1)*2+1,g[22]);    
  }
  // level 9
  c_conn[40][0] = make_comp_conn(36,g[22]);
  c_conn[40][1] = make_comp_conn(41,g[23]);
  c_conn[41][0] = make_comp_conn(36,g[17]);
  c_conn[41][1] = make_comp_conn(40,g[23]);
  c_conn[41][2] = make_comp_conn(48,g[18]);
  c_conn[41][3] = make_comp_conn(49,g[18]);
  c_conn[42][0] = make_comp_conn(37,g[17]);
  c_conn[42][1] = make_comp_conn(43,g[23]);
  c_conn[42][2] = make_comp_conn(50,g[18]);
  c_conn[42][3] = make_comp_conn(51,g[18]);
  c_conn[43][0] = make_comp_conn(37,g[22]);
  c_conn[43][1] = make_comp_conn(42,g[23]);
  c_conn[44][0] = make_comp_conn(38,g[22]);
  c_conn[44][1] = make_comp_conn(45,g[23]);
  c_conn[45][0] = make_comp_conn(38,g[17]);
  c_conn[45][1] = make_comp_conn(44,g[23]);
  c_conn[45][2] = make_comp_conn(52,g[18]);
  c_conn[45][3] = make_comp_conn(53,g[18]);
  c_conn[46][0] = make_comp_conn(39,g[17]);
  c_conn[46][1] = make_comp_conn(47,g[23]);
  c_conn[46][2] = make_comp_conn(54,g[18]);
  c_conn[46][3] = make_comp_conn(55,g[18]);
  c_conn[47][0] = make_comp_conn(39,g[22]);
  c_conn[47][1] = make_comp_conn(46,g[23]);
  for(i=0;i<4;i=i+2){
    c_conn[48+i][0] = make_comp_conn(41+i/2,g[18]);
    c_conn[48+i][1] = make_comp_conn(48+i+1,g[9]);
    c_conn[48+i][2] = make_comp_conn(56+i,g[19]);
    c_conn[48+i+1][0] = make_comp_conn(41+i/2,g[18]);
    c_conn[48+i+1][1] = make_comp_conn(48+i,g[9]);
    c_conn[48+i+1][2] = make_comp_conn(56+i+1,g[19]);
  }
  for(i=0;i<4;i=i+2){
    c_conn[52+i][0] = make_comp_conn(45+i/2,g[18]);
    c_conn[52+i][1] = make_comp_conn(52+i+1,g[9]);
    c_conn[52+i][2] = make_comp_conn(60+i,g[19]);
    c_conn[52+i+1][0] = make_comp_conn(45+i/2,g[18]);
    c_conn[52+i+1][1] = make_comp_conn(52+i,g[9]);
    c_conn[52+i+1][2] = make_comp_conn(60+i+1,g[19]);
  }
  for(i=0;i<8;i++){
    c_conn[56+i][0] = make_comp_conn(48+i,g[19]);
  }
  c_conn[64][0] = make_comp_conn(28,g[24]);
  c_conn[64][1] = make_comp_conn(65,g[25]);
  c_conn[65][0] = make_comp_conn(64,g[25]);
  c_conn[65][1] = make_comp_conn(66,g[26]);
  for(i=0;i<2;i++){
    c_conn[66+i][0] = make_comp_conn(66+i-1,g[26]);
    c_conn[66+i][1] = make_comp_conn(66+i+1,g[26]);
  }
  c_conn[68][0] = make_comp_conn(67,g[26]);
  return;
}

struct comp_conn make_comp_conn(int comp,double gamma){
  struct comp_conn c_conn;
  c_conn.comp = comp;
  c_conn.gamma = gamma;
  return c_conn;
}

double min(double a, double b){
  if(a < b){
    return a;
  }else{
    return b;
  }
}

int my_round(double x){
  if(ceil(x)-x > x-floor(x)){
    return floor(x);
  }else{
    return ceil(x);
  }
}

