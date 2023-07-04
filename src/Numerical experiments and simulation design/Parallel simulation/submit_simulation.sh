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

usage() {
    echo "Usage: submit_Excalibur.sh [OPTIONS] projectdir"
    echo ""
    echo "-h : usage"
    echo "-E : email [mandatory, must be @uclouvain.be]"
    echo "-M : MEM_job Memory per job to be used [default is 4g]"
    echo "-Q : QUEUE is the Slurm partition to use [default is cpu]"
    echo "-R : RUNNINGDIR is where to run the program, please do not change this [default is /tmp]"
    echo "-T : THREAD is the number of threads to be used per job [default is 4]"
    echo "-Y : RVERSION is the R version to be used [default is Extended-R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1]"
    echo "-Z : RLIBRARIES is the R library to be used [default is /storage/research/dduv/gehu-mv/rlib/4.2]"
    echo "-a : path_to_functions [default is /storage/home/boutrys/simulation/functions]"
    echo "-b : if you use it you will perform a power analysis, boolean to decide wether to perform a type I error OR power simulation [default is 1, which means TRUE, will perform type I error]"
    echo "-c : commonRare, boolean to decide wether to keep only rare variant or include also common variant [by default is set to 0, which means FALSE, only rare variant will be kept]"
    echo "-d : weighting_scheme, weighting procedure to apply to variant [by default is set to "Default", from seegun lee et al, can be set to "classic" or "madsen"]"
    echo "-e : causal_percent, percent of causal variant amongst variant within a region consider as causal [by default is set to 20]"
    echo "-f : negative_percent, percent of causal variant with protective effect amongst causal variant within a region [by default is set to 0]"
    echo "-g : cohort_size, size of the cohort to be analyzed, be carefull that can lead to computational explosion [by default is set to 100]"
    echo "-h : disease_preval, estimated prevalence of disease within the general population [by default is set to 0.01]"
    echo "-i : sub_region_size, genomic size of each region to be analyzed [by default is set to 3000]"
    echo "-j : prop_caseVScontrol, proportion of case within the cohort [by default is set to 0.5]"
    echo "-k : causal_maf_cutoff, MAF under which a variant might be included as causal variant [by default is set to 0.03]"
    echo "-l : Nbr_repeat, number of repeat of experiment, be carefull that my lead to computational time explosion, each repeat will be the number of jobs to be launch! recommended Nbr_repeat + 1000 for power analysis, and max up to 10 000 for type I error analysis [by default is set to 1000]"
    echo "-m : min_pos_in_sub_region, number of causal variant that has to be present in a region to be considered valid for analysis [by default is set to 2]"
    echo "-n : Type_pheno, [To Be developed] "Binary" to perform Binary analysis case versus control, "Continuous" to perform continuous analysis where phenotype is a continuous variable [by default is set to "Binary"]"
    echo "-o : Nbr_haplotype, number of haplotype to generate from the COSI input files in the litterature so far always 10 000 [by default is set to 10000]"
    echo "-p : maximum_OR, coefficient to compute weight, beta = maximum_OR * |log10MAFj| [by default is set to log(5)/4, see doi: 10.1016/j.ajhg.2011.05.029]"
    echo "-q : name_test, [To be developed] allow to run on specific aggregationt test, see testing_stat_framework.R to see name tests available [by default do not use that parameter and all test will be performed]"

    echo "-z : alpha_level, alpha level to be evaluated 0 = c(0.01, 0.05, 10^(-3)) | 1 = 0.05 | 2 = 0.01 | 3= 10^(-3) | 4 = 10^(-4) | 5 = 2.5*10^(-6) | 6 = c(0.01, 0.05, 10^(-3), 10^(-4)) | 7 = c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6)) [by default is set to 0]"

    echo "projectdir : directory of the project where to store the results. This directory must contains two tsv files named patient_data and control data (if all data for control and patient are already been retrieved) OR patient_list and control_list (if only the names of individuals should be use)((this option is mandatory)."
    echo ""
    echo "submit_simulation.sh submit the analysis of to perform simulation based on input parameters and simulated data from COSI."
    echo ""
}

## Formatting and colors (use echo -e)
RESET='\033[0m'
BOLD='\033[1m'
UL='\033[4m'
INV='\033[7m'
RED='\033[38;5;1m'
GREEN='\033[38;5;34m'
BLUE='\033[38;5;21m'
YELLOW='\033[38;5;226m'
MAGENTA='\033[38;5;165m'
PURPLE='\033[38;5;129m'
ORANGE='\033[38;5;208m'
PINK='\033[38;5;198m'
GREY='\033[38;5;246m'

if [ $# -eq 0 ]
then
    usage
    exit
fi

while getopts "hE:M:Q:R:T:Y:Z:a:bcd:e:f:g:h:i:j:k:l:m:n:o:p:q:z:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        E)
            EMAIL="$OPTARG"
            ;;
        M)
            MEM_job="$OPTARG"
            ;;
        Q)
            QUEUE="$OPTARG"
            ;;
        R)
            RUNNINGDIR="$OPTARG"
            ;;
        T)
            THREAD="$OPTARG"
            ;;
        Y)
            RVERSION="$OPTARG"
            ;;
        Z)
            RLIBRARIES="$OPTARG"
            ;;
        a)
            path_to_functions="${OPTARG}"
            ;;
        b)
            typeI=0
            ;;
        c)
            commonRare=1
            causal_maf_cutoff=0.99
            ;;
        d)
            weighting_scheme="${OPTARG}"
            ;;
        e)
            causal_percent="${OPTARG}"
            ;;
        f)
            negative_percent="${OPTARG}"
            ;;
        g)
            cohort_size="${OPTARG}"
            ;;
        h)
            disease_preval="${OPTARG}"
            ;;
        i)
            sub_region_size="${OPTARG}"
            ;;
        j)
            prop_caseVScontrol="${OPTARG}"
            ;;
        k)
            causal_maf_cutoff="${OPTARG}"
            ;;
        l)
            Nbr_repeat="${OPTARG}"
            ;;
        m)
            min_pos_in_sub_region="${OPTARG}"
            ;;
        n)
            Type_pheno="${OPTARG}"
            ;;
        o)
            Nbr_haplotype="${OPTARG}"
            ;;
        p)
            maximum_OR="${OPTARG}"
            ;;
        q)
            name_test="${OPTARG}"
            ;;
        z)
            alpha_level="${OPTARG}"
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

shift $(($OPTIND - 1))
if [ "$1" == "" ]
then
    usage
    exit
fi


#R Version and libraries
if [ "$RVERSION" == "" ]
then
    RVERSION=R-excalibur/1.0-foss-2022b-R-4.2.2
    #RVERSION=/opt/R-4.1.0-RStudio/bin/Rscript
fi
###Change here path to the right R libraries to be used
if [ "$RLIBRARIES" == "" ]
then
    RLIBRARIES=/apps/modulesfiles/software/${RVERSION}
    #RLIBRARIES="/data/gehu/rlib/4.2.0/"
fi


### store arguments
BASE_PROJECT=$(basename $1)
###Get date and time of run and add it to the name of the analysis
DATE=`date +"%Y_%m_%d_%H_%M_%S"`
PROJECT="${DATE}_${BASE_PROJECT}"
#PROJECT=${BASE_PROJECT}

### paths
SIMULATION_EXCALIBUR=/storage/home/boutrys/simulation
INPUT=${SIMULATION_EXCALIBUR}/input/${BASE_PROJECT}
LOG=${SIMULATION_EXCALIBUR}/logs/${PROJECT}
RESULTS=${SIMULATION_EXCALIBUR}/results/${PROJECT}

### We create a new folder within the logs folder for each analysis
mkdir -p ${LOG}
chmod 775 ${LOG} || true
mkdir -p ${RESULTS}
chmod 775 ${RESULTS} || true


### Retrieve error from this script
touch "${LOG}/submit_simulation.txt"
LOG_submit="${LOG}/submit_simulation.txt"
chmod 664 ${LOG_submit}
echo "assigning script variable" > ${LOG_submit}

###Default values for variables
if [ "$EMAIL" == "" ]
    then
        echo "email adress must be provided, please read usage -E"
        exit
fi
if [ "$MEM_job" == "" ]
then
    MEM_job=4g
fi
if [ "$QUEUE" == "" ]
then
    QUEUE=cpu
fi
if [ "$RUNNINGDIR" == "" ]
then
    RUNNINGDIR=/tmp
fi
if [ "$THREAD" == "" ]
then
    THREAD=4
fi
if [ "$path_to_functions" == "" ]
then
    path_to_functions="${SIMULATION_EXCALIBUR}/data/functions"
fi
if [ "$typeI" == "" ]
then
    typeI=1
fi
if [ "$commonRare" == "" ]
then
    commonRare=0
fi
if [ "$weighting_scheme" == "" ]
then
    weighting_scheme="Default"
fi
if [ "$causal_percent" == "" ]
then
    causal_percent=20
fi
if [ "$negative_percent" == "" ]
then
    negative_percent=0
fi
if [ "$cohort_size" == "" ]
then
    cohort_size=100
fi
if [ "$disease_preval" == "" ]
then
    disease_preval=0.01
fi
if [ "$sub_region_size" == "" ]
then
    sub_region_size=3000
fi
if [ "$prop_caseVScontrol" == "" ]
then
    prop_caseVScontrol=0.5
fi
if [ "$causal_maf_cutoff" == "" ]
then
    causal_maf_cutoff=0.03
fi
if [ "$Nbr_repeat" == "" ]
then
    Nbr_repeat=1000
fi
if [ "$min_pos_in_sub_region" == "" ]
then
    min_pos_in_sub_region=2
fi
if [ "$Type_pheno" == "" ]
then
    Type_pheno="Binary"
fi
if [ "$Nbr_haplotype" == "" ]
then
    Nbr_haplotype=10000
fi
if [ "$maximum_OR" == "" ]
then
    maximum_OR="classic"
fi
if [ "$name_test" == "" ]
then
    name_test=""
fi
if [ "$alpha_level" == "" ]
then
    alpha_level=0
fi

echo -e "DONE assigning script variable\n" >> ${LOG_submit}



### Submission

#maybe remove mail REQUEUE
#time=500:00:00 bigger?
#Slurm
echo -e "Submitting main job to cluster\n" >> ${LOG_submit}
sbatch --parsable --job-name=${PROJECT} --partition=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL,REQUEUE --ntasks=1 --cpus-per-task=${THREAD} --mem=${MEM_job} --time=500:00:00 --output=${LOG}/%j.%x --export=ALL,PROJECT=${PROJECT},QUEUE=${QUEUE},RUNNINGDIR=${RUNNINGDIR},EMAIL=${EMAIL},THREAD=${THREAD},MEM_job=${MEM_job},LOG=${LOG},path_to_functions=${path_to_functions},typeI=${typeI},commonRare=${commonRare},weighting_scheme=${weighting_scheme},causal_percent=${causal_percent},negative_percent=${negative_percent},cohort_size=${cohort_size},disease_preval=${disease_preval},sub_region_size=${sub_region_size},prop_caseVScontrol=${prop_caseVScontrol},causal_maf_cutoff=${causal_maf_cutoff},Nbr_repeat=${Nbr_repeat},min_pos_in_sub_region=${min_pos_in_sub_region},Type_pheno=${Type_pheno},Nbr_haplotype=${Nbr_haplotype},maximum_OR=${maximum_OR},name_test=${name_test},alpha_level=${alpha_level},LOG=${LOG},INPUT=${INPUT},RESULTS=${RESULTS},ZIP=${ZIP},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},R_LIBS=${R_LIBS} ${path_to_functions}/do_simulation.sh
echo -e "DONE Submitting main job to cluster\n" >> ${LOG_submit}


