/*
 General functions that are needed for all cells.
 */

#include "cell.h"

struct comp_conn make_comp_conn(int comp,double gamma){
  struct comp_conn c_conn;
  c_conn.comp = comp;
  c_conn.gamma = gamma;
  return c_conn;
}

struct current make_current(int comp, double (*cur_func)(int,double)){
  struct current cur;
  cur.comp = comp;
  cur.cur_func = cur_func;
  return cur;
}
