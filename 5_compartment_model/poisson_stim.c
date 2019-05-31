#include <stdlib.h>
#include <stdio.h>
#include "my_math.h"
#include "poisson_stim.h"

static int cells; // number of cells
static double *last_stim;
static double time_step;
static double stim; // amount of stimulation (typically nA)
static double stim_length;
static double lambda;
static int record_stim = 0;
static FILE* outp;

void poisson_stim_setup(int n_cells, double dt, double I, double stim_l, 
			double rate){
  extern int cells;
  extern double *last_stim;
  extern double time_step;
  extern double stim;
  extern double stim_length;

  cells = n_cells;
  last_stim = (double*)calloc(n_cells,sizeof(double));
  if(last_stim == NULL){
    fprintf(stderr,"Calloc failed to allocated memory for Poisson stim" 
	    "database.\n");
  }
  time_step = dt;
  stim = I;
  stim_length = stim_l;
  lambda = rate;
}

double poisson_stim_get_time_step(){
  extern double time_step;
  return time_step;
}

void poisson_stim_set_time_step(double dt){
  extern double time_step;
  time_step = dt;
}

double poisson_stim_get_stim(){
  extern double stim;
  return stim;
}

void poisson_stim_set_stim(double I){
  extern double stim;
  stim = I;
}

double poisson_stim_get_stim_length(){
  extern double stim_length;
  return stim_length;
}

void poisson_stim_set_stim_length(double stim_l){
  extern double stim_length;
  stim_length = stim_l;
}

double poisson_stim_get_lambda(){
  extern double lambda;
  return lambda;
}

void poisson_stim_start_recorder(char *filename, char *mode){
  extern int record_stim;
  extern FILE* outp;

  record_stim = 1;
  outp = fopen(filename,mode);
}

void poisson_stim_stop_recorder(){
  extern int record_stim;
  extern FILE* outp;

  record_stim = 0;
  fclose(outp);
}

void poisson_stim_reset(){
  // reset internal database to 0
  extern double *last_stim;
  int i;

  for(i=0;i<cells;i++){
    last_stim[i] = 0;
  }
}

void poisson_stim_dismantle_database(){
  extern double *last_stim;
  free(last_stim);
}

double poisson_stim(int id, double t){
  extern double *last_stim;
  extern double stim;
  extern FILE* outp;
  extern int record_stim;
  extern double stim_length;
  
  if((double)rand()/RAND_MAX < time_step*lambda*exp(-time_step*lambda)){
    if(record_stim){
      fprintf(outp,"%d\t%g\n",id,t);
    }
    *(last_stim+id) = t;
    return stim;
  }else if(last_stim[id] != 0 && t-last_stim[id] <= stim_length){
    return stim;
  }else{
    return 0;
  }
}
