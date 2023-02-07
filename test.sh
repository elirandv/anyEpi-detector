ANYEPI_ENV_NAME="5gmc_n3_detector"
echo "activate ${ANYEPI_ENV_NAME} conda environment"
echo SHELL_NAME
conda init bash
conda activate 5gmc_n3_detector #$ANYEPI_ENV_NAME

echo "activate ${ANYEPI_ENV_NAME} conda environment"
conda activate $ANYEPI_ENV_NAME
if [$CONDA_DEFAULT_ENV!=$ANYEPI_ENV_NAME]; then
  raise error "Error: ${ANYEPI_ENV_NAME} conda environment activation Failed!"
echo "activate ${ANYEPI_ENV_NAME}"

exit