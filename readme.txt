This is the readme for the models associated with the paper:

Munro E, Borgers C (2010) Mechanisms of very fast oscillations in
networks of axons coupled by gap junctions. J Comput Neurosci

These model files were contributed by Erin Munro
ecmun at bu.edu

Example build and run: Under linux an executable can be built in the
5_compartment_model folder with the following command:

gcc -lm Traub1999_5comp.c CA3pyramidal_axon.c cell.c gap_junction.c \
        my_math.c connections.c poisson_stim.c -o run.exe

For code runnability a test can be started with a command like
./run.exe -V -65 -60 -gj 10 10.1 -seed 1 10 -run_time 1 -stim_stop 0.5
(this is only for demonstration of the format of the command)

Details about the model files and parameters are provided below.

5_compartment_model:

This folder contains all the code for the 5-compartment axon
model. The model is written in C, but uses an object-oriented
style. The axon model is contained in CA3pyramidal_axon.c, which is
the 5 axonal compartments from Traub et al. (1999) with a fixed
somatic voltage. The axon uses the framework from cell.c, where you
can make compartment connections and add currents. The two currents
that can be added are gap junctional currents (gap_junction.c) and
external Poisson stimulation (poisson_stim.c).

Traub1999_5comp.c simulates the large network model of an axonal
plexus using the axon from CA3pyramidal_axon.c and the same network as
Traub et al. (1999) (implemented in connections.c). It takes the
following arguments:

-V V_Sl V_Su: run simulation for somatic voltages V_Sl to V_Su (0.2 mV
between each voltage)
-gj g_gjl g_gju: run simulation for gap junction conductances g_gjl to
g_gju (0.1 nS between each conductance)
-seed seedl seedu: run simulation for seeds seedl to seedu (seeds 1-10
used in Munro & Borgers (2010)
-run_time: length of simulation in ms
-stim_stop: time to stop Poisson stimulation in ms

For each (V_s, g_gj) pair, it prints out a data file that is readable
into MATLAB as a matrix. The first column gives the time and each
subsequent column gives the voltage (in mV) of an axon. Note, not all
calculated times are listed in this matrix for efficiency.

graph_refractory_test.c simulates the small network model and outputs
the time in between the two stimulations (delta) along with the firing
times of axons 0, 1, and 2 (1-3 in Munro & Borgers (2010)). It takes
the following arguments:
-V V_Sl V_Su: run simulation for somatic voltages V_Sl to V_Su (0.2 mV
between each voltage)
-gj g_gjl g_gju: run simulation for gap junction conductances g_gjl to
g_gju (0.1 nS between each conductance)
-delta deltal deltau: run simulation for delta values deltal to deltau
(dt=0.0025 ms between each delta)

Traub1999_reproduction:
This folder contains our implementation in C of the model for figure
12 in Traub et al. (1999). The CA3 pyramidal cell model is in
traub69.c. The network is implemented in connections.c. run_traub.c
puts these together to run the full simulation. The output is printed
in data files readable by MATLAB as matrices. Each data file contains
500 ms of data. The first column stands for time. Columns 2 and 3
stand for the voltages of the soma and second-to-last axonal
compartment of the cell 1. Likewise, all subsequent columns list the
voltages of the soma and then second-to-last axonal compartment of
each cell in order. Note, not all calculated times are printed in
these files for efficiency reasons.

cellular_automata:
This folder contains the MATLAB codes for the two modified cellular
automata to model noise and re-entrant activity from the axonal
plexus. The file conn.out contains the network generated for the
models in C (all axonal plexus simulations and the Traub 1999
reproduction use this network).

CA_sol_double.m contains the general cellular automaton to model
noise, and run_CA_double_exp.m runs the experiment to mimick noise.

CA_sol_ref.m contains the general cellular automaton to model
re-entrant activity, and run_CA_ref_exp.m runs the experiment to
mimick re-entrant activity with various refractory periods for
4-connected cells.

MATLAB_files:

Here are some MATLAB files that were used to analyze the data,
including analysis on large random networks as well as processing of
data from the 5-compartment model.

cconn2mconn.m - converts connections network generated in C to cell
array readable by other MATLAB files
cluster_placement.m - arranges cells according to cluster size
conn_analysis.m - prints out brief synopsis of properties of network,
and returns the clusters within the network
count_4_conn_neighbors.m - counts down stream neighbors of 4-connected
cells to determine "typical" gap junctional load
find_cycles.m - gives an underestimate of the number of cycles going
through a cell
freq.m - generates the power spectrum and plots the results
remove_cells.m - returns a modified network where all connections to
4-connected cells are removed
spike_failures.m - lists propagations failures (and successes) in a
given data set
spike_timing.m - lists the spike times from a given data set and plots
a rastergram
