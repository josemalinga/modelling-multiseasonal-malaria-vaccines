Basic workflow for creating scenario files and running OpenMalaria simulations

created by Monica Golumbeanu and modified by lydia burgert
monica.golumbeanu@unibas.ch, lydia.burgert@unibas.ch
12.02.2021
__________________________________

This workflow is built to run from https://rstudio.scicore.unibas.ch/. 

DESCRIPTION:
The workflow started with the OM_base_workflow_original.R can be used to specify the paramters for new OM simulations, create scenario .xml files and run OpenMalaria simulations for each scenario.

Functions:
- create_folders.R: creates the necessary folders in user and GROUP M3TPP folder 
-generate_param_table.R: generetates parameter table by using lhs sampling to draw continuous parameter samples and then combining with all combinations of categorical variables

-genOMsimscripts: generates the scripts job_create_scenarios.sh, run_OM.sh and submission_*expname*_*no*.sh scripts
The main script of the workflow is  submission_*expname*_*no*.sh and calls sequentially two array job scripts: submit_create_scenario.sh which creates all the scenarios and run_OM.sh which runs OpenMalaria simulations for all the scenarios, resource requirements and job diagnostics need to be specified in this script by replacing the corresponding options in the job specification header
#!/bin/bash
#SBATCH --job-name=run_RCD
#SBATCH --account=penny
#SBATCH -o /scicore/home/penny/*user*/M3TPP/*exp*/JOB_OUT/*scriptname*.out
#SBATCH --mem=1G
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1

REQUIRED RESOURCES:
OpenMalaria resource files needed for running the workflow are located in the folder OM_schema38/

You can change the OpenMalaria version by modifying the load module command in the script run_OM.sh
In addition, you will also need to copy the corresponding files as provided in OM_schema38/ for the new version
You can specify a new location for the resource folder in the script run_OM.sh

INPUT:
To run the workflow, you need to 
0) create an "Experiments" folder within the M3TPP project in your home directory
1) make a folder in M3TPP/Experiments with your experiment name. Then copy the file OM_base_workflow_original.R in analysisworkflow into that experiment folder and give it the suffix of your experiment name (replace original to your name).
2)Inside this script specify the experiment name and chunk size (batch size for simulation runs, can not be over 100 000) and run create_folders.R (creates folders OM_JOBS and JOB_OUT).
3)Then copy a scaffold.xml file into the OM_JOBS folder which contains marked parameters (surrounded by @) which will be varied across simulations. Check 0_scenarios/ for resources. File needs to be re-named scaffold.xml
4) copy a seasonality.txt file specifying the seasonality profiles into the ./M3TPP/Experiments/"exp"/ folder. Check 0_scenarios/ for resources. File needs to be renamed seasonlity.txt
5)specify the parameters you want to sample (continious) and categorical parameters for scenario definition. Parameter names need to correspond to the parameters marked with @ in the xml file.
6) specify number of samples and number of seeds and run the functions gen_paramtable.R and genOMsimscripts.R.



OUTPUTS:
The workflow creates the following sub-folders located in the generated GROUP folder:

base/ contains for each scenario replacement patterns for the parameters 
scenarios/ contains the final scenario.xml files to be used for openMalaria simulations
om/ contains the results of the OpenMalaria simulations




