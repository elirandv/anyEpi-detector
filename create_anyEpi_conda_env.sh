ANYEPI_ENV_NAME=$1

CUDA_VERSION=$(nvidia-smis | grep -oP "CUDA Version: \K.*"  | grep -oE "\-?[0-9]+\.[0-9]+")
echo "CUDA Version on nvidia-smi:${CUDA_VERSION}"

if [ -z "$CUDA_VERSION" ]
then
  echo "no CUDA installed, please install CUDA Version>11.0"
  exit
elif [ $CUDA_VERSION < 11 ]
then
  echo "Bad CUDA Version, please install CUDA Version>11.0"
  exit
fi

echo "create ${ANYEPI_ENV_NAME} conda environment with environment.yml dependencies"
conda env create --name $ANYEPI_ENV_NAME --file environment.yml

echo "activate ${ANYEPI_ENV_NAME} conda environment"
conda activate $ANYEPI_ENV_NAME

HDF5_PLUGIN_PATH=$(whereis hdf5)
export HDF5_PLUGIN_PATH=$HDF5_PLUGIN_PATH/lib/plugin/

#echo "install requirements.txt"
#pip install -r requirements.txt
#
#
#echo "add env paths"
#
#echo "GPU test results: ${}"