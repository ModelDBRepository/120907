#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include "connections.h"

#define PI 3.14159265358979

void traub69(int,int[][4],double,int,double,int*,int,int);

int main(int argc, char **argv){
  double dt = .0025;
  double run_time = 100;
  int start_rec = 0;
  int n_row = 32;
  int n_col = 96;
  int cells = n_row*n_col;
  int nRanks,rank,n_per_rank,n_first,n_last,leftover,cell_assgn[81],
    conn[cells][4];
  double r_c = 10;
  double c = 1.6;
  int i,k;
  FILE* outp;
  struct rlimit rlim;
  struct rusage usage;
  int error;
  double V_S,g_gj,stim,lambda;
  int seed;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("argc = %d\n",argc);
  for(i=0;i<argc;i++){
    printf("argv[%d] = %s\n",i,argv[i]);
  }
  fflush(NULL);

  n_per_rank = cells/nRanks;
  leftover = fmod(cells,nRanks);
  n_first = 0;
  for(i=0;i<nRanks;i++){
    cell_assgn[i] = n_first;
    if(leftover>0){
      leftover = leftover - 1;
      n_first = n_first + n_per_rank + 1;
    }else{
      n_first = n_first + n_per_rank;
    }
  }
  cell_assgn[nRanks] = cells;

  connections(conn,n_row,n_col,r_c,c);
  
  if(rank == 0){
    outp = fopen("conn.out","w");
    for(i=0;i<(n_row*n_col);i++){
      for(k=0;k<4;k++){
        fprintf(outp,"%d\t",conn[i][k]);
      }
      fprintf(outp,"\n");
    }  
    fclose(outp);
  }
  
  getrlimit(RLIMIT_STACK,&rlim);
  
  rlim.rlim_cur = rlim.rlim_max;
  setrlimit(RLIMIT_STACK,&rlim);
  
  //defaults
  g_gj = 3.66;
  stim = .2;
  lambda = .002; // stimulus/ms/cell %1
  seed = 0;
	    
  traub69(cells,conn,run_time,0,dt,cell_assgn,rank,nRanks)

  MPI_Finalize();
  return 0;
}

