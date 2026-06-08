#!/bin/bash


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
mydir="$(cd `dirname $0` && pwd)"
source "${mydir}/bash_install_utils.sh"

# Check that we don't have conflicting environmental variables set.
if pytool_cfg_conflict; then
    cat << MSG
Both GGG_USE_PIP and GGG_USE_MICROMAMBA are set to 1 or greater.
Please ensure only one of these variables is set to 1 or greater
to select an install tool.
MSG
    exit 2
fi


# Try pip first, although this isn't the preferred approach. As long as
# we leave conda for last as the default, this script will work since we
# checked above that there aren't conflicted environmental variables.
if do_use_pip_only; then
    PYTHON_CMD="$GGG_PYTHON"
    if [ -z "$PYTHON_CMD" ]; then
        PYTHON_CMD=$(which python3)
    fi
    if [ -z "$PYTHON_CMD" ]; then
        cat << MSG

python3 not found, ensure python3 is on your path or set the GGG_PYTHON
environmental variable to the full path to your python program.
MSG
        exit 2
    fi
    if [ ! -e "$PYTHON_CMD" ]; then
        cat << MSG

Python does not exist at path "$PYTHON_CMD".
If you have the GGG_PYTHON environmental variable set, ensure it points to
an extant python program. Otherwise, please check the output of "which python3"
to determine why it is pointing to a nonexistant file.
MSG
        exit 2
    fi
    pyversion=$("$PYTHON_CMD" --version | awk '{print $2}')
    if version_ge $pyversion "3.10.0"; then
        echo "Found good Python version $pyversion"
    else
        cat << MSG

The Python version given (${pyversion}) is less than 3.10.0.
GGG2020.1 only supports Python 3.10 and up. Please provide the
path to a more recent version of python as the GGG_PYTHON
environmental variable.
MSG
        exit 2
    fi


    echo "export ENVDIR='$GGGPATH/install/.venv'" > $mydir/.init_conda
    echo "export NO_CONDA_ACTIVATE=1" >> $mydir/.init_conda
    echo "export PY_FOR_CREATE='$PYTHON_CMD'" >> $mydir/.init_conda
    echo "export PIPCMD='$mydir/.venv/bin/pip'" >> $mydir/.init_conda
    echo "export PYTHONCMD='$mydir/.venv/bin/python'" >> $mydir/.init_conda
    exit 0
fi

# Try micromamba next - since it's a standalone executable, we shouldn't
# need all of the shell init stuff like for conda below
if do_use_micromamba; then
  which micromamba > /dev/null
  mm_found=$?
  if [ $mm_found == 0 ] || [ $GGG_USE_MICROMAMBA == 2 ]; then
    echo "export ENVDIR='$GGGPATH/install/.condaenv'" > $mydir/.init_conda
    echo "CONDACMD='micromamba --yes'" >> $mydir/.init_conda
    echo "CREATECMD='env create'" >> $mydir/.init_conda
    echo "UPDATECMD='update'" >> $mydir/.init_conda
    echo "export PYTHONCMD='micromamba run -p $mydir/.condaenv python'" >> $mydir/.init_conda
    echo "export PIPCMD='micromamba run -p $mydir/.condaenv pip'" >> $mydir/.init_conda
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

https://tccon-wiki.caltech.edu/Main/InstallingPython

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
    echo "export ENVDIR='$GGGPATH/install/.condaenv'" >> $mydir/.init_conda
    echo "CONDACMD=conda" >> $mydir/.init_conda
    echo "CREATECMD='env create'" >> $mydir/.init_conda
    echo "UPDATECMD='env update'" >> $mydir/.init_conda
    echo "export PYTHONCMD=python" >> $mydir/.init_conda
    echo "export PIPCMD=pip" >> $mydir/.init_conda
    echo "export NO_CONDA_ACTIVATE=0" >> $mydir/.init_conda
fi

