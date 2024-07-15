#!/bin/bash -l

always_yes=false

for arg in $@; do
    case $arg in 
        -y|--yes)
            always_yes=true
            ;;
    esac
done

# Can't use $GGGPATH because the login shell (-l in the shebang)
# will restore GGGPATH to default, so if we've changed it, we 
# have problems
mydir=$(cd `dirname $0` && pwd)

# Try micromamba first - since it's a standalone executable, we shouldn't
# need all of the shell init stuff like for conda below
if [ ! -z $GGG_USE_MICROMAMBA ] && [ $GGG_USE_MICROMAMBA -gt 0 ]; then
  which micromamba > /dev/null
  mm_found=$?
  if [ $mm_found == 0 ] || [ $GGG_USE_MICROMAMBA == 2 ]; then
    echo "CONDACMD='micromamba --yes'" > $mydir/.init_conda
    echo "CREATECMD='env create'" >> $mydir/.init_conda
    echo "UPDATECMD='update'" >> $mydir/.init_conda
    echo "export PYTHONCMD='micromamba run -p $mydir/.condaenv python'" >> $mydir/.init_conda
    echo "export NO_CONDA_ACTIVATE=1" >> $mydir/.init_conda
    exit 0
  else
    echo "micromamba not on path, trying conda"
  fi
fi
  

conda_base="$(conda info --base 2>/dev/null)" || conda_base=false
if [ $conda_base == false ] ; then
        cat << MSG 

Anaconda3 is not installed
To install Anaconda3, you may use the install_conda.sh script or 
manually install/download. If you already have Anaconda installed 
but it is not being detected, visit

https://tccon-wiki.caltech.edu/Main/RunningGGG2020#Python_47Anaconda

for advice to debug this problem.

MSG
        exit 2
elif $always_yes ; then
        echo "Will use anaconda3 in $conda_base"
else
        echo "Anaconda3 exists in: $conda_base"
        read -p "Use the Anaconda installed in $conda_base ? [yn]" answer
        case $answer in
                [yY])
                        ;;
                *)
                        cat <<MSG

You have chosen not to use the Anaconda installed at 

$conda_base

An Anaconda install is required for GGG. If you have a different
install you wish to use, make sure that 'conda' called from a bash
login shell resolves to the desired installation.

MSG
                        exit 2
                        ;;
        esac
fi

# Build a hidden file in the same directory as this file to source to activate conda
conda_init_file=${conda_base}/etc/profile.d/conda.sh
if [ ! -f $conda_init_file ]; then
    echo "Cannot find conda init file at: $conda_init_file"
    exit 1
else
    echo "source $conda_init_file" > $mydir/.init_conda
    echo "CONDACMD=conda" >> $mydir/.init_conda
    echo "CREATECMD='env create'" >> $mydir/.init_conda
    echo "UPDATECMD='env update'" >> $mydir/.init_conda
    echo "export PYTHONCMD=python" >> $mydir/.init_conda
    echo "export NO_CONDA_ACTIVATE=0" >> $mydir/.init_conda
fi

