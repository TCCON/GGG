#!/bin/bash

# Necessary to use conda commands inside this script, even though this was already 
# done in master.sh. This will set CONDACMD to either "conda" or "micromamba"
source "$GGGPATH/install/.init_conda"


if [ -e "$GGGPATH/install/.condaenv" ]; then
        echo "The $GGGPATH/install/.condaenv python environment already exists; update from environment.yml"
        $CONDACMD $UPDATECMD -p "$GGGPATH/install/.condaenv" -f environment.yml
else
        echo "The $GGGPATH/install/.condaenv python environment does not exist"
        echo "Creating python environment from environment.yml"
        $CONDACMD $CREATECMD -p "$GGGPATH/install/.condaenv" -f environment.yml
fi
