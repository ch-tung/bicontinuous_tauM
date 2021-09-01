* `in.maxwell_long`: `LAMMPS` input file.  Generates the stress tensor as a function of time (`maxwell_long.out`) .
* `load_thermo.py`: Reads `maxwell_long.out` and transform the data into `.mat` format `thermo*.mat`.
* `tau_M_long.py`: Calculates the stress autocorrelation (`maxwell_long_py.mat`).
* `relaxation_time.m`: Loads `maxwell_long_py.mat` and calculated the Maxwell relaxation time. Generate figures of stress autocorrelation and tau_M versus temperature.

