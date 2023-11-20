# Shell_Mass
Code for Shell Mass Estimation paper "Galactic HI supershells: kinetic energies and possible origin", LA Suad, CF Caiafa, EM Arnal, S Cichowolski, Astronomy&Astrophysics, Vol. 624, Apr 2019. doi: https://doi.org/10.1051/0004-6361/201833850

To generate the results reported in the paper the scripts named "Compute_Xc_Dvel_YYY.m", located at the /Scripts_paper/, folder need to be run, where X is the quadrant (3 or 4) and YYY is the Delta Velocity.

The information of previously identified shells are included in the functions "load_shells_3c_Dvel_100.m", "load_shells_3c_total.m", "load_shells_4c_Dvel_100.m" and "load_shells_4c_total.m".

Some relevant functions are included in the folder /main_functions/ such as:

- "run_all_masses.m": This function compute the masses and temperatures of each of the shells.
- "Compute_mass_V5.m": This function computes the radial profiles of a shell starting from its central local minima. Based on these profiles, the external walls are detected. The cloud of points detected as potential walls are then refined and finally an ellipse is fit to them.


