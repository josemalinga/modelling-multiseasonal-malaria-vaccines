# modelling-multiseasonal-malaria-vaccines

Performance characteristics and potential public health impact of improved malaria vaccines targeting childhood burden

Josephine Malinga1,2, Lydia Braunack-Mayer1,2, Thiery Masserey1,2, Matthew Cairns3, Sherrie L Kelly1,2, Narimane Nekkab1,2, Melissa A Penny1,2

1 Swiss Tropical and Public Health Institute, Allschwil, Switzerland 2 University of Basel, Basel, Switzerland 3 London School of Hygiene and Tropical Medicine, London, United Kingdom 

*Correspondence to: melissa.penny@unibas.ch

Using an individual-based malaria transmission model (https://github.com/SwissTPH/openmalaria/wiki) we evaluate the population level impact of deploying a long duration pre_erythrocytic vaccine, linking vaccine performance properties, health system and programmatic factors and understanding their trad-offs and relationships.

Folders / Workflow Steps
analysisworkflow
This workflow builds on the workflow presented in Golumbeanu (2021) and Burgert (2021) to specify Target Product Profiles for new interventions against malaria. First, a set of simulated scenarios is defined. These are characterized by the delivery modality, tool specifications, and settings in which a concrete health target is analysed. Second, a set of disease scenarios are simulated randomly over the entire parameter space to evaluate the health outcomes. The resulting database of simulations is used to train a Gaussian process emulator (GP), that predicts the health outcome given a set of input parameters. Third, the emulator is employed to perform sensitivity analysis and optimisation of tool properties with respect to health outcomes. This analysis allows to define the optimal product characteristics of new interventions that maximises the chance of achieving a desired health goal.

Contributors (in chronological order): Melissa Penny, Guojing Yang, Monica Golumbeanu, Lydia Burgert, Mirjam Laager, Narimane Nekkab, Josephine Malinga, Lydia Braunack-Mayer

0_scenarios
Contains XML files and associated parameter ranges used to simulate data with OpenMalaria (https://github.com/SwissTPH/openmalaria/wiki).
1_OM_basic_workflow
Generates paramater table and XML scenarios from base scaffold.xml
Launches OM simulations with 2 outputs
The ctsout.txt contains a table with outputs for "continuous time" the measures specified in the the scenario test.xml. There is one line for each (5-day) time step.
The output.txt contains a table with four columns and no headers for survey measures.
2_postprocessing
Performs generalized post-processing of OM simulations by the settings specified in previous sets
For each setting, a split file is generated in “postprocessing/split” that specifies the parameters for this setting and based on that, a seeds file (for every simulation) and an agg file (aggregated over all seeds for one parameter set) is generated
For each iTPP, postprocessing functions have been further developed and saved in 0_scenarios
3_GP_train
Trains GPs using for a specified outcome calculated in 2_postprocessing for each of the seeds files.
predictors: continuous variables
predicted: health outcome
To change GP plotting and outputs modify script analysisworkflow/3_GP_train/train_GP.R. Adaptive sampling can be added in this step depending on GP performance.
4_sensitivity_analysis
Performs sensitivity analysis of health outcome using analysis of variance (sobol) for GPs trained in step 3 within the parameter bounds used for simulation (default)
predictors: continuous variables
predicted: health outcome
To chance number of bootstrap samples, change function calc_sobol_idx in analysisworkflow/3_GP_train/GPtoolbox.R
5_optimization
Performs non-linear optimisation of chosen continuous input variables to reach a certain health goal while keeping other continuous variables constant
If the grid size is to wide, change number of grid points within 5_optimization/optimize_parameter.R
The non-linear search method performs optimisation by using the Augmented Lagrange Multiplier Method with a pre-trained emulator in “3_GP_train”.
6_grid_optimization
Alterative to step 5, performs a grid search optimisation of chosen continuous input variables to reach a certain health goal while keeping other continuous variables constant
The grid search method uses a pre-trained emulator in “3_GP_train” for optimisation.
data_and_visualisation
This folder contains the data generated during this study, along with the R scripts used to visualise data. There is a folder for each figure in the manuscript and supplement, containing:

The .rds data file corresponding to the figure,
The Rscript used to generate the figure, and
A pdf version of the figure.
To reproduce a given figure, download the corresponding folder and update the file paths referenced in the corresponding Rscript.
