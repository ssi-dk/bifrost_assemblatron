#!/bin/bash
helpFunction()
{
   echo ""
   echo "Usage:"
   echo "$0 -i LOCAL - for local install"
   echo "$0 -i COMP - for computerome install"
   exit 1 # Exit script after printing help
}
while getopts "i:" opt
do
   case "$opt" in
      i ) parameterI="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
# Print helpFunction in case parameters are empty
if [ -z "$parameterI" ]
then
   echo "-i parameter is empty";
   helpFunction
fi
if [ "$parameterI" != "LOCAL" ] && [ "$parameterI" != "COMP" ]
then
  echo "Wrong argument"
  helpFunction
fi
if [ "$parameterI" == "LOCAL" ]
then
  echo "Starting local install"
fi
if [ "$parameterI" == "COMP" ]
then
  echo "Starting computerome install"
  module load tools computerome_utils/2.0
  module load tools anaconda3/2022.10
fi
# Begin script
if $(conda config --show channels | grep -q "bioconda")
then
  echo "bioconda channel is already added"
else
  conda config --add channels bioconda
  echo "bioconda channel was added, you can remove it after installation with command"
  echo "conda config --remove channels bioconda"
fi
if $(conda config --show channels | grep -q "conda-forge")
then
  echo "conda-forge channel is already added"
else
  conda config --add channels conda-forge
  echo "bioconda channel was added, you can remove it after installation with command"
  echo "conda config --remove channels conda-forge"
fi
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REQ_TXT="$SCRIPT_DIR/environment.yml"
CONFIG_YAML_PATH=$(find $SCRIPT_DIR -name "config.yaml")
if test -f "$CONFIG_YAML_PATH";
then
  COMPONENT_NAME=$(grep "display_name:.*." $CONFIG_YAML_PATH | tr " " "\n" | grep -v "display_name:")
  if [ -z "$COMPONENT_NAME" ]
    then
      echo "display_name: in config.yaml should contain component name"
      exit 1
  fi
  COMPONENT_VERSION=$(grep -o "code:.*." $CONFIG_YAML_PATH | tr " " "\n" | grep -v "code:")
  if [ -z "$COMPONENT_VERSION" ]
    then
      echo "code: in config.yaml should contain component version"
      exit 1
  fi
  ENV_NAME=("bifrost_"$COMPONENT_NAME"_"$COMPONENT_VERSION)
else
  echo "Cannot find config.yaml in component folder to form env name"
  exit 1
fi
if test -f "$REQ_TXT";
then
  echo "Making conda env"
  echo "$ENV_NAME will be created"
  conda env create --file "$REQ_TXT" --name $ENV_NAME
else
  echo "environment.yml file cannot be found in the script folder"
  exit 1
fi
