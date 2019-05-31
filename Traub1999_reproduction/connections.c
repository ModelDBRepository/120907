/* Establishes connections for n_row by n_col grid with
   max Euclidean distance r_c, and c connections
   per cell placed randomly. conn is an (n_row*n_col,1)
   array. Row i gives the numbers of the cells cell i
   is connected to.
*/
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "connections.h"

int within_bounds(int*,int*,int,int,double);
int already_connected(int conn[][4],int,int);

void connections(int conn[][4], int n_row,int n_col,double r_c, double c){
  
  int x,y,i,a[2],b[2],a_index,b_index,rel_coord;
  unsigned long m, max_conn;
  int M=n_row*n_col;
  int n_conn[M];
  int d = floor(r_c);
  int rc_list[(int)pow(2*d+1,2)][2];

  // initialize the conn and n_conn matrices
  for(x=0;x<M;x++){
    n_conn[x] = 0;
    for(i=0;i<4;i++){
      conn[x][i] = -1;
    }
  }
  
  srand(2);
  
  // number of connections to make
  if(M>1){
    max_conn = factorial(M-1);
  }else{
    max_conn = 0;
  }
  if(ceil(M*c/2) < max_conn){
    m = ceil(M*c/2);
  }else{
    m = max_conn;
  }
  
  // relative indexes of cells r_c from (0,0)
  i = 0;
  for(x=-d;x<=d;x++){
    for(y=-d;y<=d;y++){
      rc_list[i][0] = x;
      rc_list[i][1] = y;
      i = i+1;
    }
  }

  i = 0;
  while(m>0){
    // choose random cell a
    a_index = floor((double)rand()*M/RAND_MAX);
    index_to_coord(a_index,n_row,n_col,a);
    
    // choose random cell b, at most r_c away from a,
    // include phantom cells
    rel_coord = floor(rand()*pow(2*d+1,2)/RAND_MAX);
    b[0] = rc_list[rel_coord][0] + a[0];
    b[1] = rc_list[rel_coord][1] + a[1];
    b_index = coord_to_index(b[0],b[1],n_row,n_col);
    
    if(within_bounds(a,b,n_row,n_col,r_c) & 
       n_conn[a_index] < 4 & n_conn[b_index] < 4 & 
       !already_connected(conn,a_index,b_index)){
      // cell is not phantom, connection not already made
      // set connection
      conn[a_index][n_conn[a_index]] = b_index;
      n_conn[a_index] = n_conn[a_index] + 1;
      conn[b_index][n_conn[b_index]] = a_index;
      n_conn[b_index] = n_conn[b_index] + 1;;
      m = m-1;
    }
    i = i+1;
  }
  return;
}

int coord_to_index(int x,int y,int n_row,int n_col){
  // translates row x, col y to cell #,
  // ordered by row then column
  // first row/col is #0, first cell is #0
  
  int index = x*n_col + y;
  return index;
}

int* index_to_coord(int index,int n_row,int n_col,int* coord){
  // give the coordinates in row and column # when
  // given the cells #
  // first row/col is #0, first cell is #0
  
  // row coordinate
  coord[0] = index/n_col;

  // column coordinate
  coord[1] = fmod(index,n_col);
  return coord;
}

int within_bounds(int* a,int* b,int n_row,int n_col,double r_c){
  /* Check if b is on the grid and in range of a */
  int yes;
  yes = 0<=b[0] & b[0]<n_row;
  yes = yes & 0<=b[1] & b[1]<n_col;
  yes = yes & sqrt(pow(a[0]-b[0],2)+pow(a[1]-b[1],2))<=r_c;
  return yes;
}

int already_connected(int conn[][4], int a_index,int b_index){
  int yes = 0;
  int k;
  if(a_index == b_index){
    return 1;
  }
  for(k=0;k<4;k++){
    if(conn[a_index][k] == b_index){
      yes = 1;
    }
  }
  return yes;
}

unsigned long factorial(int N){
  int  n;
  unsigned long N_fact;
  n = N;
  N_fact = 1;
  while(n>1){
    if(N_fact > ULONG_MAX/n){
      return ULONG_MAX;
    }
    N_fact = N_fact*n;
    n = n-1;
  }

  return N_fact;
}
