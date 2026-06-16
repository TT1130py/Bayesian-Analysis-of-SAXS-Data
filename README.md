# Bayesian Analysis of SAXS Data
### Ensemble Refinement of Intrinsically Disordered Proteins against small angle x-ray scattering data

***

#### Intrinsically disordered proteins are a class of biomolecules that exist as a highly dynamic ensemble of different conformations. Computational methods only provide static representations of IDP structures, while experimental techniques can only record ensemble-averaged information. To examine the three dimensionanal structure dynamics of an IDP ensemble, a Bayesian Inference derived reweighting protocol to reweight probable structures based on their fit to experimental data.

***
### Path recommendations
#### Most of this code depends on specific file directories linked to certain files and paths that contain results from individual parts of the main iBME simulation. The directories can be modified to the user's liking; however, for ease of use, the following directory tree should be followed to implement with the default code.

```
my-project/  
├── analysis/  
|   └── curve_compare.py                     #all of the analysis code in this folder 
├── grid_sims                                #this will be the main output folder for all grid sim runs before reweighting
|   ├── saxs_dro_x_to_y_r0_x2_to_y2/         #individual grid sim run defined by the drho and r0 values
|       ├── compiled_GPs                     #compiled saxs curves for every grid point from iBME_reweight
|       ├── GPx                              #data for each grid point
|       ├── iBME_results                     #output from iBME_reweight.py- this folder should eventually be moved to its own directory
|       └── grid_full.txt                    #the full grid data
├── iBME_res/                                #main output folder for all iBME reweight runs
|   ├── iBME_dro_x_to_y_r0_x2_to_y2_theta_z/ #individual reweight run defined by drho, r0, and theta values
|       ├── GRID_sum.txt                     #summary of fitting parameters for every gp
|       └── any other output files that can be ran from the analysis foler should go here
├── BME.py                                   #main iBME framework
├── BME_tools.py                             #functions for running iBME
├── iBME_script.py                           #edited function script for running the iBME method
├── do_gp_fraction.sh                        #shell script needed in the main pipeline
├── Pepsi-SAXS                               #actual Pepsi SAXS executable
├── iBME_env/                                #venv
└── all of the main pipeline scripts
```
***






