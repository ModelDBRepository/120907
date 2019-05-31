void gap_junction_setup(int n_cells, int (*connections)[4],
			 double conductance);
double* gap_junction_database();
double gap_junction_conductance();
void gap_junction_set_conductance(double conductance);
void gap_junction_dismantle_database();
double gap_junction_current(int id, double t);
