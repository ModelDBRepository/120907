/*
 General functions that are needed for all cells.
 */

/* a compartment connection meant to go in an array indexed by primary 
   compartments then the comp_conn is a compartment connected to the primary 
   compartment with conductance gamma */
struct comp_conn {
  int comp;
  double gamma;
};

/* a current that is added to the right hand side of the differential 
   equation for compartment comp. The current is given by cur_func which 
   takes a cell id and a time (t) */
struct current {
  int comp;
  double (*cur_func)(int,double);
};

struct comp_conn make_comp_conn(int comp,double gamma);
struct current make_current(int comp, double (*cur_func)(int,double));
