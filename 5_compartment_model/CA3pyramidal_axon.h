struct CA3pyramidal_axon {
  int id;
  // variables for all compartments
  double V_S;
  double V[5];
  double m[5];
  double h[5];
  double n[5];
  struct current *currents; // pointer to array of currents
  int n_currents;
};

void CA3pyraxon_setup_constants();
void CA3pyraxon_init(struct CA3pyramidal_axon *axon, int id, double V_S);
void CA3pyraxon_set_currents(struct CA3pyramidal_axon* axon, 
			     struct current *currents, 
			     int n_currents);
void CA3pyraxon_step(struct CA3pyramidal_axon *axon0, 
		   struct CA3pyramidal_axon *axon1,
		   double dt, double t);
