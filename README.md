# How to run everything
1. Run A1/scripts/aero_design.m

Saves variables to aero_design.mat

2. Run A1/scripts/HAWC_inputs.m
Creates operation\_\*.dat files for:
  - Single point operation
  - Multiple TSR operation
  - Multiple wind speed operation

Creates data/redesign/redesign_ae.dat

Creates A1/data/redesign_c2def.dat

Prints genspeed for operational_data block

3. Run \*\_aero_rigid.htc with HAWC2S

4. Move results to A1/res/\*/tsr or tsr_opt or ws

5. Run A1/scripts/postProcessing_aero.m

6. Run A2/scripts/Structural_scaling.m
  - Writes data/redesign/redesign_Blade_st.dat
  - Writes data/redesign/redesign_c2def.dat
  - Saves variables to mat/Structural_scaling.mat

Now, you either:
  - Create info/body.dat and info/beam.dat

Uncomment lines 5 and 6 and running HAWC2S

(*NB: Change line 103 nbodies = 1 for blades),*

(NB: you don't need an operational_data_filename for this)

  - Create .opt file

Comment out info/body.dat and info/beam.dat

Comment in compute_optimal_pitch_angle

Run with HAWC2S

  - Create results

Leave info/body.dat and info/beam.dat commented out

Comment in operational_data_filename

Comment in compute_steady_state

Comment in whatever save_\* you want

Move results to A2/res/\*

7. Run A2/scripts/postProcessing_struct.m

Save variables to mat/postProcessing_struct.mat

8. Create Campbell and damping diagrams with HAWCStab2

Instructions in Section 5.3 of wind_turbine_modal_analysis.pdf

Save Campbell diagrams (\*\_ael.cmb, \*\_st.cmb) to A2/res/\*/stab/

Save binary results (\*\_ael.hsd, \*\_st.hmd) to A2/res/\*/stab/

where \* is turbine (DTU_10MW or redesign)

9. Run A2/scripts/postProcessing_stab.m

Writes mat/postProcessing_stab.mat

10. Run redesign_cont_HS2.htc using HAWC2S

Check genspeeds (*NB: gearratio = 1 here)*

Change full load omega and zeta in controller_tuning block

Choose constant power or constant torque and gain scheduling fit

10. Run A3/scripts/postProcessing_cont_abridged.m

Writes data/redesign/controller_gains.dat

*11. Check USB/hawc_simulations/htc_master/redesign_DL.htc*

e.g. minimum rotor speed, scaling of aerodynamic output channels (spanwise position), box_dim_u and box_dim_v

12. Run turbulent simulations

See notes/hawc_simulations.txt

13. Move res and log from USB to corresponding directory in A4/

14. Run A4/scripts/post_process_hawc2.py for turbulent simulations

15. Evaluate performance using Python scripts

- calculate_aep.py (saves comparison bar chart to A4/figs/Part4/AEP.png)
- calculate_extreme_loads.py (saves polar plot to A4/figs/Part4/polar_extr.png)
- calculate_fatigue_loads.py (saves polar plot to A4/figs/Part4/polar_fat.png)
- plot_dels_2turbines.py (saves short-term fatigue DELs to A4/figs/Part3/)
- plot_maxmeanmin_2turbines.py (saves max, mean, min plots to A4/figs/Part1and2/)

*NB: remember to change the wind class. I.e. For AEP both turbines should be in class IIIb, for loads each turbine should be in the wind class it was designed for, i.e. DTU 10MW in Ia and redesign in IIIb*

NB: DTU 10MW results folder split into turbulence class a (tca/) and turbulence class b (tcb/) results
