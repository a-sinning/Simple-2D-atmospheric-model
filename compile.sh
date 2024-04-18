rm -f allroutines.f90 onefile.f90 model.exe *.o *.mod *out.txt 

cat driver.f90 setup_outfiles.f90 terrain_setup.f90 terrainBCs.f90 terrain_filt_horiz.f90 terrain_filt_vert.f90 physical_grid_locs.f90 read_sounding.f90 compute_basestate.f90 writeout_basestate.f90 initial_rand_noise.f90 flux_predictive_eqns.f90 spatial_filtering.f90 saturation_adjust.f90 modeltop_damping.f90 timefilter.f90 handoff_values.f90 zerogradLB_pert_var_BCs.f90 periodicLB_pert_var_BCs.f90 leapfrog_setup.f90 warmbubble.f90 writeout_perts.f90 close_files.f90 > allroutines.f90

gfortran -fbounds-check -o model.exe usersettings.f90 declarations.f90 allroutines.f90

cat usersettings.f90 declarations.f90 allroutines.f90 > onefile.f90
