#!/bin/bash
#
#  Set shell, otherwise the default shell would be used.
#$ -S /bin/bash
#
#  Set memory limit.
#$ -l h_vmem=10G
#
#  Set wallclock time limit - depending on batch size.
#$ -l h_rt=3:00:00
#
#  Make sure that the .e and .o file arrive in the
#  working directory.
#$ -cwd
#
#  Specify path to job's stdout output file.
#$ -o ../logs/Cityscapes/trainval_refined_stereogfSemanticRefinement_beta_0.01
#
#  Merge the standard out and standard error to one file.
#$ -j y
#
#  Declare job name
#$ -N Fog_Cityscapes
#
#  Force / switch off an immediate execution of the job:
#  crucially affects execution success
#$ -now n
#
#  Schedule K jobs with ids 1-K. The specific value of K
#  depends on the batch size.
#$ -t 1-11
#
source /home/sgeadmin/BIWICELL/common/settings.sh
/bin/echo Running on host: `hostname`
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`
/bin/echo PATH: `echo $PATH`
/bin/echo TMP: `env | grep TMP` 
#/bin/echo SGE: `env | grep SGE`
/bin/echo MCR: `env | grep MCR`

# Parameters.
TASK_ID=${SGE_TASK_ID:-"$1"}
DATASET_SPLIT="trainval"
REFINEMENT_LEVEL="refined"
VARIANT="stereoscopic_inpainting_with_guided_filtering_and_semantic_segment_refinement"
BETA="0.01"
OUTPUT_ROOT_DIR="/scratch_net/nowin/csakarid/Toyota-foggy/data/Cityscapes" # Change to preferred directory for writing results of simulation.
IMAGES_PER_TASK="50"

# Change directory to that containing the fog simulation script.
cd ..

# Fog simulation script.
/usr/sepp/bin/matlab -nodesktop -nodisplay -nosplash -r "Fog_simulation_Cityscapes('${TASK_ID}', '${DATASET_SPLIT}', '${REFINEMENT_LEVEL}', '${VARIANT}', ${BETA}, '${OUTPUT_ROOT_DIR}', ${IMAGES_PER_TASK}); exit;"

echo finished at: `date`
exit 0;
