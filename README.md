To equilibrate the system for low density (in this case \phi=0.01) is necessary the follow procedure:
  1. Generate the cubic lattice for macroions such that the box is small enough to accomodate the conterions. Run lattice_cubic.c code.
  2. At this time the counterions are on spherical layer, to do this configuration we run sphere.c code. In this code change the number of counterions around the macroion.
  3. After, run the Monte Carlo simulation to equilibrate the macroion and counterions system.
  4. To produce dynamical properties run the Brownian Dynamic simulations.
