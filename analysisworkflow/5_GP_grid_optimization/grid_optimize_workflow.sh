#SBATCH --qos=1day

SIM_FOLDER=$1
PRED=$2
SCALE=$3
NGRID=$4
TARGET_RANGE_SIZE=$5

PARAM_RANGES_FILE=$SIM_FOLDER"param_ranges.RData"
GP_FOLDER=$SIM_FOLDER"gp/trained/"$PRED"/"

# echo "SIM_FOLDER: $SIM_FOLDER"
# echo "PRED: $PRED"
# echo "SCALE: $SCALE"
# echo "NGRID: $NGRID"
# echo "TARGET_RANGE_SIZE: $TARGET_RANGE_SIZE"

# Submit GP optimization analysis array job
gp_files=(${GP_FOLDER}*.RData)
NUM=${#gp_files[@]}

sbatch --array=1-$NUM job_grid_optimize_parameter.sh $GP_FOLDER $SIM_FOLDER $PRED $SCALE $NGRID $TARGET_RANGE_SIZE

