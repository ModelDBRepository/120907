/*
  Using a network with one center axon with 4 neighbors, stimulates the first 
  neighbor twice for delta=0:dt:10 to graph the refractory period of the 
  center cell. The refractory period is the shortest length of time after the 
  center axon fires before it can fire again. (arrays notated as in MATLAB)
*/

#include "gap_junction.h"
#include <stdio.h>
#include <stdlib.h>
#include "my_math.h"
#include "cell.h"
#include "CA3pyramidal_axon.h"

#define n_cells 14 // total number of cells
#define rec_cells 3 // number of cells to be recorded

static struct CA3pyramidal_axon cells[2][n_cells];
static double delta;
static double results[rec_cells][2];
// for each cell - time of first firing, time of second firing,
// 0 stands for never firing
static int up[rec_cells];
static struct current crnt[2];
static double run_time = 30;
static double dt = .0025;

void run(double V_S);
double test_stim(int id, double t);

int main(int argc, char* argv[]){
  extern struct CA3pyramidal_axon cells[2][n_cells];
  extern double delta;
  extern struct current crnt[2];
  extern double dt;
  int conn[n_cells][4];
  double V_Sl,V_Su; // mV
  double g_gjl,g_gju; // nS
  double deltal,deltau; // ms
  double V_S = 0; // mV
  double g_gj = 4; // nS
  int i,j,k,m,n_conn;
  int n[3];
  char file[100];
  FILE* outp;
  
  i = 0;
  while(i<argc){
    if(0==strcmp("-n",argv[i])){
      // number of connections for cell 9
      // should be 3 or 4
      n[0] = atoi(argv[i+1]);
      // number of connections per cell: 10-12
      // can be 1-4
      n[1] = atoi(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-V",argv[i])){
      V_Sl = atof(argv[i+1]);
      V_Su = atof(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-gj",argv[i])){
      g_gjl = atof(argv[i+1]);
      g_gju = atof(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-delta",argv[i])){
      deltal = atof(argv[i+1]);
      deltau = atof(argv[i+2]);
      i = i+3;
    }else{
      i = i+1;
    }
  }   

  for(i=0;i<n_cells;i++){
    for(j=0;j<4;j++){
      conn[i][j] = -1;
    }
  }
  conn[0][0] = 1;
  conn[1][0] = 0;
  if(n[0]>=3){
    conn[1][1] = 2;
    conn[1][2] = 3;
    conn[2][0] = 1;
    conn[3][0] = 1;
    if(n[0]==4){
	conn[1][3] = 4;
	conn[4][0] = 1;
    }
  }
  if(n[1]>=2){
    for(i=0;i<3;i++){
      conn[2+i][1] = 2+3*(i+1);
      conn[2+3*(i+1)][0] = 2+i;
    }
    if(n[1]>=3){
      for(i=0;i<3;i++){
	conn[2+i][2] = 2+3*(i+1)+1;
	conn[2+3*(i+1)+1][0] = 2+i;
      }
      if(n[1]>=4){
	for(i=0;i<3;i++){
	  conn[2+i][3] = 2+3*(i+1)+2;
	  conn[2+3*(i+1)+2][0] = 2+i;
	}
      }
    }
  }

  //printf("created network\n");
  //fflush(NULL);
  
  CA3pyraxon_setup_constants();
  crnt[0] = make_current(4,(double (*)(int,double))test_stim);
  gap_junction_setup(n_cells,conn,g_gj*1e-3);
  crnt[1] = make_current(3,(double (*)(int,double))gap_junction_current);
  for(i=0;i<2;i++){
    CA3pyraxon_init(&cells[i][0],0,V_S);
    CA3pyraxon_set_currents(&cells[i][0],crnt,2);
    for(k=1;k<n_cells;k++){
      CA3pyraxon_init(&cells[i][k],k,0);
      CA3pyraxon_set_currents(&cells[i][k],&crnt[1],1);
    }
  }
  
  
  for(i=0;i<=my_round((V_Su-V_Sl)/.2);i++){
    V_S = V_Sl+.2*i;
    printf("V_S = %g\n",V_S);
    for(j=0;j<=my_round((g_gju-g_gjl)/.1);j++){
      g_gj = g_gjl+.1*j;
      printf("g_gj = %g\n",g_gj);
      gap_junction_set_conductance(g_gj*1e-3);

      //sprintf(file,"/cluster/shared/emunro01b/grf_extended/grf_extended_VS%dgj%dn%d-%d.out",
      sprintf(file,"grf_extended_VS%dgj%dn%d-%d.out",
	      my_round(V_S*10),my_round(gap_junction_conductance()*1e4),n[0],n[1]);
      outp = fopen(file,"w");
      printf("starting runs\n");
      fflush(NULL);
      //for(i=0;i<=ceil(10/dt);i++){
      for(k=0;k<=my_round((deltau-deltal)/dt);k++){//ceil(10/dt);k++){
	delta = deltal + dt*k;//k*dt;
	//printf("delta = %g\n",delta);
	//fflush(NULL);
	run(V_S);
	fprintf(outp,"%g",delta);
	for(m=0;m<rec_cells;m++){
	  fprintf(outp,"\t%g\t%g",results[m][0],results[m][1]);
	}
	fprintf(outp,"\n");
	fflush(NULL);
      }
      fclose(outp);
    }
  }

      
  gap_junction_dismantle_database();
}

void run(double V_S){
  extern struct CA3pyramidal_axon cells[2][n_cells];
  double *pac;
  extern double results[rec_cells][2];
  extern int up[rec_cells];
  extern struct current crnt[2];
  extern double run_time;
  extern double dt;
  struct CA3pyramidal_axon *old_pyr, *new_pyr;
  int i,j;

  pac = gap_junction_database();

  for(i=0;i<rec_cells;i++){
    up[i] = 0;
    for(j=0;j<2;j++){
      results[i][j] = 0;
    }
  }

  for(i=0;i<2;i++){
    for(j=0;j<n_cells;j++){
      CA3pyraxon_init(&cells[i][j],j,V_S);
    }
    CA3pyraxon_set_currents(&cells[i][0],crnt,2);
    for(j=1;j<n_cells;j++){
      CA3pyraxon_set_currents(&cells[i][j],&crnt[1],1);
    }
  }

  for(i=1;i<=ceil(run_time/dt);i++){

    // swap new cells to old position
    if(fmod(i,2)==0){
      old_pyr = &cells[1][0];
      new_pyr = &cells[0][0];
    }else{
      old_pyr = &cells[0][0];
      new_pyr = &cells[1][0];
    }

    // update all cells via midpoint method
    for(j=0;j<n_cells;j++){
      CA3pyraxon_step(&old_pyr[j],&new_pyr[j],dt,i*dt);
    }

    for(j=0;j<n_cells;j++){
      // update voltages for gj's all at once for consistency
      pac[j] = new_pyr[j].V[3];
    }

    for(j=0;j<rec_cells;j++){
      if(up[j]==0 && new_pyr[j].V[3] > 50){
        // just starting to spike
        if(results[j][0] == 0){
          results[j][0] = i*dt;
        }else if(results[j][1] == 0){
          results[j][1] = i*dt;
          if(j==rec_cells-1){
            return;
          }
        }
        up[j] = 1;
      }else if(up[j]==1 && new_pyr[j].V[3] <= 50){
        up[j] = 0;
      }
    }

  }
  return;
}

double test_stim(int id, double t){
  extern double delta;

  if(id==0 && 10<t && t<10.3125){
    return 0.2;
  }else if(id==0 && 10+delta<t && t<10.3125+delta){
    return 0.2;
  }else{
    return 0;
  }
}

