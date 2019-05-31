#include <stdlib.h>
#include <stdio.h>
#include "connections.h"
#include "gap_junction.h"

static double g_gj; // gap junctions conductance in uS
static double *pac; // pointer to array of voltage of penultimate axonal compartment of everyone
static int (*conn)[4]; // gap junction connections (max 4 conn, -1 if no conn)

void gap_junction_setup(int n_cells, int (*connections)[4],double conductance){
  extern double g_gj;
  extern double *pac;
  int i,j;

  g_gj = conductance;
  pac = (double*)calloc(n_cells,sizeof(double));
  if(pac == NULL){
    fprintf(stderr,"Calloc failed to allocated memory for gap junctions"
            "database.\n");
  }
  conn = connections;
}

double* gap_junction_database(){
  // needed to update voltages
  extern double *pac;
  return pac;
}

double gap_junction_conductance(){
  extern double g_gj;
  return g_gj;
}

void gap_junction_set_conductance(double conductance){
  extern double g_gj;
  g_gj = conductance;
}

void gap_junction_dismantle_database(){
  extern double *pac;
  free(pac);
}

double gap_junction_current(int id, double t){
  extern double *pac;
  extern int (*conn)[4];
  extern double g_gj;
  double I = 0;
  int j = 0;
  
  while(j<4 && conn[id][j] != -1){
    //printf("Adding current from %d to %d\n",conn[id][j],id);
    I = I + g_gj*(pac[conn[id][j]] - pac[id]);
    j = j+1;
  }
  
  return I;
}
