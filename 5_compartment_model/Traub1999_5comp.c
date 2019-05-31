/*
  The same simulation from Traub 1999, only using the 5 compartments of the
  axon and initial segment with a fixed somatic voltage.
*/

#include "gap_junction.h"
#include <stdio.h>
#include <stdlib.h>
#include "my_math.h"
#include "cell.h"
#include "connections.h"
#include "CA3pyramidal_axon.h"
#include "poisson_stim.h"

#define n_row 32
#define n_coll 96
#define n_cells 3072
static struct CA3pyramidal_axon cells[2][n_cells];
static struct current crnt[2];
static double run_time = 100;
static double stim_stop_time = 50;
static double dt = .0025;
static double V_S;
static double g_gj;
static int seed = 1;

void run(double V_S);

int main(int argc, char* argv[]){
  extern struct CA3pyramidal_axon cells[2][n_cells];
  extern struct current crnt[2];
  extern double dt;
  int conn[n_cells][4];
  double V_Sl,V_Su; // mV
  double g_gjl,g_gju; // nS
  extern double V_S,g_gj;
  extern double run_time;
  extern double stim_stop_time;
  int i,k,j,seedl,seedu;
  extern int seed;
  
  i = 0;
  while(i<argc){
    if(0==strcmp("-V",argv[i])){
      V_Sl = atof(argv[i+1]);
      V_Su = atof(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-gj",argv[i])){
      g_gjl = atof(argv[i+1]);
      g_gju = atof(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-seed",argv[i])){
      seedl = atoi(argv[i+1]);
      seedu = atoi(argv[i+2]);
      i = i+3;
    }else if(0==strcmp("-run_time",argv[i])){
      run_time = atof(argv[i+1]);
      i = i+2;
    }else if(0==strcmp("-stim_stop",argv[i])){
      stim_stop_time = atof(argv[i+1]);
      i = i+2;
    }else{
      i = i+1;
    }
  }   
  
  connections(conn,32,96,10,1.6);
  
  CA3pyraxon_setup_constants();
  // must use dt/2 since using midpoint method
  poisson_stim_setup(n_cells,dt/2,.2,.3125,.002);
  crnt[0] = make_current(4,(double (*)(int,double))poisson_stim);
  gap_junction_setup(n_cells,conn,g_gjl*1e-3);
  crnt[1] = make_current(3,(double (*)(int,double))gap_junction_current);
  for(i=0;i<2;i++){
    printf("crnt[%d] = {%d, %p}\n",i,crnt[i].comp,crnt[i].cur_func);
  }

  for(i=0;i<=my_round((V_Su-V_Sl)/.2);i++){
    V_S = V_Sl+.2*i;
    printf("V_S = %g\n",V_S);
    for(j=0;j<=my_round((g_gju-g_gjl)/.1);j++){
      g_gj = g_gjl+.1*j;
      printf("g_gj = %g\n",g_gj);
      gap_junction_set_conductance(g_gj*1e-3);
      for(k=0;k<=(seedu-seedl);k++){
	seed = seedl+k;
	printf("seed = %d\n",seed);
	poisson_stim_reset();
	srand(seed);
	run(V_S);
      }
    }
  }
  
  poisson_stim_dismantle_database();
  gap_junction_dismantle_database();
}

void run(double V_S){
  extern struct CA3pyramidal_axon cells[2][n_cells];
  double *pac;
  extern double results[n_cells][2];
  extern struct current crnt[2];
  extern double run_time;
  extern double dt;
  struct CA3pyramidal_axon *old_pyr, *new_pyr;
  double rec_step;
  int i,j,k;
  char file[100];
  FILE* outp;

  //sprintf(file,"/scratch/emunro01b/Traub1999_5comp_VS%dgj%dseed%d.out",
  sprintf(file,"Traub1999_5comp_VS%dgj%dseed%d.out",
	  my_round(V_S*10),my_round(g_gj*10),seed);
  //sprintf(file,"Traub1999_5comp_VS%dgj%dseed%d.out",
  //  my_round(V_S*10),my_round(g_gj*10),seed);
  outp = fopen(file,"w");
  rec_step = my_round(.1/dt);
  
  pac = gap_junction_database();

  for(i=0;i<2;i++){
    for(j=0;j<n_cells;j++){
      CA3pyraxon_init(&cells[i][j],j,V_S);
      CA3pyraxon_set_currents(&cells[i][j],crnt,2);
    }
  }

  printf("stim_stop_time = %g, ceil(stim_stop_time/dt) = %g\n",
	 stim_stop_time, ceil(stim_stop_time/dt));
  fflush(NULL);
  sprintf(file,"/scratch/emunro01/Traub1999_5comp_VS%dgj%dseed%d_stim.out",
	  my_round(V_S*10),my_round(g_gj*10),seed);
  //poisson_stim_start_recorder(file,"w");
  for(i=1;i<=ceil(stim_stop_time/dt);i++){

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
    
    if(fmod(i,rec_step)==0){
      fprintf(outp,"%g\t",i*dt);
      for(j=0;j<n_cells;j++){
	fprintf(outp,"%g\t",new_pyr[j].V[3]);
	fprintf(outp,"%g\t",new_pyr[j].m[3]);
	fprintf(outp,"%g\t",new_pyr[j].h[3]);
	fprintf(outp,"%g\t",new_pyr[j].n[3]);
      }
      fprintf(outp,"\n");
    }
  }
  //poisson_stim_stop_recorder();

  for(i=0;i<2;i++){
    for(j=0;j<n_cells;j++){
      CA3pyraxon_set_currents(&cells[i][j],&crnt[1],1);
    }
  }

  printf("run_time = %g, ceil(run_time/dt) = %g\n",run_time,
	 ceil(run_time/dt));
  for(i=ceil(stim_stop_time/dt)+1;i<=ceil(run_time/dt);i++){

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
    
    if(fmod(i,rec_step)==0){
      fprintf(outp,"%g\t",i*dt);
      for(j=0;j<n_cells;j++){
	fprintf(outp,"%g\t",new_pyr[j].V[3]);
	fprintf(outp,"%g\t",new_pyr[j].m[3]);
	fprintf(outp,"%g\t",new_pyr[j].h[3]);
	fprintf(outp,"%g\t",new_pyr[j].n[3]);
      }
      fprintf(outp,"\n");
    }
  }

  fclose(outp);
  return;
}
