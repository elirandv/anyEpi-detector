EPICALL_ENV_NAME="EpiCall"
CONDA_ENVS_DIR="/ebensteinlab/eliraneitan/condaenvs/"

#CUDA_VERSION=$(nvidia-smis | grep -oP "CUDA Version: \K.*"  | grep -oE "\-?[0-9]+\.[0-9]+")
#echo "CUDA Version on nvidia-smi:${CUDA_VERSION}"
#
#if [ -z "$CUDA_VERSION" ]
#then
#  echo "no CUDA installed, please install CUDA Version>11.0"
#  exit
#elif [ $CUDA_VERSION < 11 ]
#then
#  echo "Bad CUDA Version, please install CUDA Version>11.0"
#  exit
#fi

echo "module load miniconda/miniconda3-4.7.12-environmentally:"
module load miniconda/miniconda3-4.7.12-environmentally

echo "conda env create ${EPICALL_ENV_NAME} conda environment with environment.yml dependencies:"
conda env create --name $EPICALL_ENV_NAME --prefix ${CONDA_ENVS_DIR} --file environment.yml

echo "conda activate ${CONDA_ENVS_DIR}${EPICALL_ENV_NAME} conda environment:"
conda activate $CONDA_ENVS_DIR$EPICALL_ENV_NAME

#echo "install requirements.txt"
#pip install -r requirements.txt
#
#
#echo "add env paths"
#
#echo "GPU test results: ${}"