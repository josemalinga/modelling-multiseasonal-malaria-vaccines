# UPDATE IN PROGRESS

# Performance characteristics and potential public health impact of improved pre-erythrocytic malaria vaccines targeting childhood burden

**Josephine Malinga<sup>1,2</sup>, Lydia Braunack-Mayer<sup>3,4</sup>, Thiery Masserey<sup>3,4</sup>, Matthew Cairns<sup>3</sup> , Sherrie L Kelly<sup>1,2</sup>, Epke Le Rutte>1</sup>, Narimane Nekkab<sup>3,4</sup>, Melissa A Penny<sup>1,2</sup>**

_<sup>1</sup> Telethon Kids Institute, Nedlands, WA, Australia; <sup>2</sup> Centre for Child Health Research, University of Western Australia, Crawley, WA, Australia; <sup>3</sup>Swiss Tropical and Public Health Institute, Allschwil, Switzerland; <sup>4</sup>University of Basel, Basel, Switzerland; <sup>5</sup>London School of Hygiene and Tropical Medicine, London, United Kingdom_ 

*Correspondence to: melissa.penny@telethonkids.org.au

Using an individual-based malaria transmission model (https://github.com/SwissTPH/openmalaria/wiki) we evaluate the population level impact of deploying a long duration pre-erythrocytic vaccine, linking vaccine performance properties, health system and programmatic factors and understanding their trade-offs and relationships.

# Workflow Steps
## analysisworkflow
This workflow builds on the workflow presented in Golumbeanu (2021) and Burgert (2021) to specify Target Product Profiles for new interventions against malaria. First, a set of simulated scenarios is defined. These are characterized by the delivery modality, tool specifications, and settings in which a concrete health target is analysed. Second, a set of disease scenarios are simulated randomly over the entire parameter space to evaluate the health outcomes. The resulting database of simulations is used to train a Gaussian process emulator (GP), that predicts the health outcome given a set of input parameters. Third, the emulator is employed to perform sensitivity analysis and optimisation of tool properties with respect to health outcomes. This analysis allows to define the optimal product characteristics of new interventions that maximises the chance of achieving a desired health goal.

Contributors (in chronological order): Melissa Penny, Guojing Yang, Monica Golumbeanu, Lydia Burgert, Mirjam Laager, Narimane Nekkab, Josephine Malinga, Lydia Braunack-Mayer

### 0_scenarios
Contains XML files and associated parameter ranges used to simulate data with OpenMalaria (https://github.com/SwissTPH/openmalaria/wiki).

### 1_OM_basic_workflow
Generates paramater table and XML scenarios from base scaffold.xml
Launches OM simulations with 2 outputs
* The ctsout.txt contains a table with outputs for "continuous time" the measures specified in the the scenario test.xml. There is one line for each (5-day) time step.
* The output.txt contains a table with four columns and no headers for survey measures.

### 2_postprocessing
Performs generalized post-processing of OM simulations by the settings specified in previous sets
For each setting, a split file is generated in “postprocessing/split” that specifies the parameters for this setting and based on that, a seeds file (for every simulation) and an agg file (aggregated over all seeds for one parameter set) is generated
For each iTPP, postprocessing functions have been further developed and saved in 0_scenarios

### 3_GP_train
Trains GPs using for a specified outcome calculated in 2_postprocessing for each of the seeds files.
predictors: continuous variables
predicted: health outcome
To change GP plotting and outputs modify script analysisworkflow/3_GP_train/train_GP.R. Adaptive sampling can be added in this step depending on GP performance.

### 4_sensitivity_analysis
Performs sensitivity analysis of health outcome using analysis of variance (sobol) for GPs trained in step 3 within the parameter bounds used for simulation (default)
predictors: continuous variables
predicted: health outcome

To chance number of bootstrap samples, change function calc_sobol_idx in analysisworkflow/3_GP_train/GPtoolbox.R

### 5_optimization
Performs non-linear optimisation of chosen continuous input variables to reach a certain health goal while keeping other continuous variables constant
If the grid size is to wide, change number of grid points within 5_optimization/optimize_parameter.R
The non-linear search method performs optimisation by using the Augmented Lagrange Multiplier Method with a pre-trained emulator in “3_GP_train”.

### 6_grid_optimization
Alterative to step 5, performs a grid search optimisation of chosen continuous input variables to reach a certain health goal while keeping other continuous variables constant
The grid search method uses a pre-trained emulator in “3_GP_train” for optimisation.

## data_and_figures
This folder contains the data generated during this study, along with the R scripts used to visualise data. There is a folder for each figure in the manuscript and supplement, containing: 
* The .rds data file corresponding to the figure,
* The Rscript used to generate the figure, and
* A pdf version of the figure. \n
To reproduce a given figure, download the corresponding folder and update the file paths referenced in the corresponding Rscript.

# In-silico modelling to validate booster efficacy
## booster_validation
Contains a simplified workflow(as shown above), datasets and associated results files used in the in-silico modelling exercise to validate RTS,S parameters in a seasonal use case

*  ### 0_scenarios (see above)
*  ### 1_OM_basic_workflow (see above)
*  ### 2_postprocessing_validation
*  ### 3_results_and_figures