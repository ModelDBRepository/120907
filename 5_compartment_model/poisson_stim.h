void poisson_stim_setup(int n_cells, double dt, double I, double stim_l, 
			double lambda);
double poisson_stim_get_time_step();
void poisson_stim_set_time_step(double dt);
double poisson_stim_get_stim();
void poisson_stim_set_stim(double I);
double poisson_stim_get_stim_length();
void poisson_stim_set_stim_length(double stim_l);
double poisson_stim_get_lambda();
void poisson_stim_start_recorder(char *filename, char *mode);
void poisson_stim_stop_recorder();
void poisson_stim_reset();
double poisson_stim(int id, double t);
