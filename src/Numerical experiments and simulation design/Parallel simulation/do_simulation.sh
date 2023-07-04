#!/bin/bash
#/*****************************************************************************************
#  *
#  * Excalibur simulation- Copyright (C) <2017-2023> <Université catholique de Louvain (UCLouvain)>
#  * 	
#  * List of the contributors to the development of Excalibur simulation: see LICENSE file.
#* Description and complete License: see LICENSE file.
#* 	
#  * This program (Excalibur simulation) is free software: 
#  * you can redistribute it and/or modify it under the terms of the 
#* GNU General Public License as published by the Free Software Foundation, 
#* either version 3 of the License, or (at your option) any later version.
#* 
#  * This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#* 
#  * You should have received a copy of the GNU General Public License
#* along with this program (see COPYING file).  If not, 
#* see <http://www.gnu.org/licenses/>.
#* 
#  *****************************************************************************************/
#  
#  /**
#  *
#  * @author Simon Boutry
#*
#  */
###Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
#set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized


### The following mandatory variables are set in submit_Excalibur.sh and are passed as environment variables by Slurm
#FOR SLURM JOBS
###PROJECT
###QUEUE
###RUNNINGDIR
###EMAIL
###THREAD
###MEM_job
###INPUT
###RESULTS
###LOG
###ZIP
##RVERSION
##RLIBRARIES
###R_LIBS

#FOR R SCRIPT
###path_to_functions
###typeI
###commonRare
###weighting_scheme
###causal_percent
###negative_percent
###cohort_size
###disease_preval
###sub_region_size
###prop_caseVScontrol
###causal_maf_cutoff
###Nbr_repeat
###min_pos_in_sub_region
###Type_pheno
###Nbr_haplotype
###maximum_OR
###name_test
###alpha_level


### Functions 
## Output starting time, host name and stage name (stages are big chunk of the pipeline)
# First argument is a string describing the stage
function startStage()
{
	host=$(hostname)
	stageStart=$(date +%s)
	echo "### Started on $host at $(date +"%Y-%m-%d %H:%M") - $1"
}

## Output ending time, elpased time, host name and stage name (stages are big chunk of the pipeline)
# First argument is a string describing the stage (use the same than startStage)
function finishStage()
{
    host=$(hostname)
	stageEnd=$(date +%s)
	stageDiff=$(((stageEnd-stageStart)/60))
	echo "### Ended on $host at $(date +"%Y-%m-%d %H:%M") - $1"
	echo "### Elapsed on $host : $stageDiff minutes - $1"
}

## Load modules given after first argument and output starting time, host name and step name (steps are generally one command)
# First argument is a string describing the step
# You can give any number of modules to load starting as second argument
function startStep()
{
	host=$(hostname)
	stepStart=$(date +%s)
	echo "## Started on $host at $(date +"%Y-%m-%d %H:%M") - $1"
    if [ "$#" -gt "1" ]
    then
        module load ${@:2}	
    fi	
}

## Purge modules and output ending time, elpased time, host name and step name (steps are generally one command)
# First argument is a string describing the step (use the same than startStep)
function finishStep()
{
	module purge
	host=$(hostname)
	stepEnd=$(date +%s)
	stepDiff=$(((stepEnd-stepStart)/60))
	echo "## Ended on $host at $(date +"%Y-%m-%d %H:%M") - $1"
	echo "## Elapsed on $host : $stepDiff minutes - $1"
}


### Retrieve error from this script
touch "${LOG}/do_simulation.txt"
LOG_script="${LOG}/do_simulation.txt"
chmod 664 ${LOG_script}


### Module version
echo "Start of the script : Loading R module" > ${LOG_script}
export R_LIBS=/apps/modulesfiles/software/${RVERSION}
echo -e "DONE Loading R module\n" >> ${LOG_script}


### Create temporary folders within the results folder of the analysis 
echo -e "Creating temporary directories within results folder\n" >> ${LOG_script}
RESULTS_TMP="${RESULTS}/tmp"
mkdir -p ${RESULTS_TMP}
chmod 777 ${RESULTS_TMP} || true

RESULTS_TMP_INPUT="${RESULTS_TMP}/input_data"
mkdir -p ${RESULTS_TMP_INPUT}
chmod 777 ${RESULTS_TMP_INPUT} || true

RESULTS_TMP_OUTPUT="${RESULTS_TMP}/output_data"
mkdir -p ${RESULTS_TMP_OUTPUT}
chmod 777 ${RESULTS_TMP_OUTPUT} || true
echo -e "DONE Creating temporary directories within results folder\n" >> ${LOG_script}


### Data Preparation
echo -e "Data Preparation\n" >> ${LOG_script}
startStep "Data Preparation" $RVERSION
R_SCRIPT=${path_to_functions}/data_preparation.R
Rscript ${R_SCRIPT} ${path_to_functions} ${INPUT} ${RESULTS_TMP_INPUT} ${typeI} ${commonRare} ${weighting_scheme} ${causal_percent} ${negative_percent} ${cohort_size} ${disease_preval} ${sub_region_size} ${prop_caseVScontrol} ${causal_maf_cutoff} ${Nbr_repeat} ${min_pos_in_sub_region} ${Type_pheno} ${Nbr_haplotype} ${maximum_OR} ${RLIBRARIES} ${name_test}
finishStep "Data Preparation"
echo -e "DONE Data Preparation\n" >> ${LOG_script}


### Submission
echo -e "Submission of jobs\n" >> ${LOG_script}
SUB_LOG="${LOG}/region_logs"
mkdir -p ${SUB_LOG}
chmod 777 ${SUB_LOG} || true

startStage "Submitting jobs for a region analysis"
startStep "Region analysis"
cd ${RESULTS_TMP_INPUT}
nbr_region=(*)
nbr_region=${#nbr_region[@]}
files=$(basename "${RESULTS_TMP_INPUT}/*")
finishStep "Region analysis"
LIST_JOB=""
for i in ${files}
do
    namefile="${i}"
    IFS="_"
    read -ra id_region <<< "${namefile}"
    IFS=""
    #Variable for Slurm
    single_region="${PROJECT}_Region_nbr_${id_region}"
    startStep "launching ${single_region}"
    JOB_single_region=`sbatch --parsable --job-name=${single_region} --partition=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=REQUEUE --ntasks=1 --cpus-per-task=${THREAD} --mem=${MEM_job} --time=500:00:00 --output=${SUB_LOG}/%j.%x --export=ALL,id_region=${id_region},RESULTS_TMP_INPUT=${RESULTS_TMP_INPUT},RESULTS_TMP_OUTPUT=${RESULTS_TMP_OUTPUT},path_to_functions=${path_to_functions},Type_pheno=${Type_pheno},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},name_test=${name_test} ${path_to_functions}/single_region_analysis.sh`
    LIST_JOB="${LIST_JOB}:${JOB_single_region}"
    finishStep "launching ${single_region}"
    #sleep 1
done
finishStage "Submitting jobs for a region analysis"
echo -e "DONE Submission of jobs\n" >> ${LOG_script}


### Results combination
echo -e "Results combination\n" >> ${LOG_script}
#Variable for Slurm
DEPEND="afterok${LIST_JOB}"
result_combi="${PROJECT}_result_combi"
startStep "launching ${result_combi}"
JOB_result_combi=`sbatch --parsable --job-name=${result_combi} --partition=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=END,FAIL,REQUEUE --ntasks=1 --cpus-per-task=${THREAD} --mem=${MEM_job} --time=500:00:00 --output=${LOG}/%j.%x --dependency=${DEPEND} --export=ALL,RESULTS_TMP=${RESULTS_TMP},RESULTS=${RESULTS},alpha_level=${alpha_level},Nbr_repeat=${Nbr_repeat},Type_pheno=${Type_pheno},typeI=${typeI},RESULTS_TMP_OUTPUT=${RESULTS_TMP_OUTPUT},path_to_functions=${path_to_functions},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},name_test=${name_test} ${path_to_functions}/results_combination.sh`
finishStep "launching ${result_combi}"
echo -e "DONE Results combination\n" >> ${LOG_script}

