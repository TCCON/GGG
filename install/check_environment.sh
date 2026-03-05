#!/bin/bash

# Necessary to use conda commands inside this script, even though this was already 
# done in run_pyinstall.sh. This will set CONDACMD to either "conda" or "micromamba"
source "$GGGPATH/install/.init_conda"

install_needed=y
md5_file="$GGGPATH/install/.condaenv/environment.md5"
if [[ -e "$GGGPATH/install/.condaenv" ]]; then
        echo "The $GGGPATH/install/.condaenv python environment already exists; verifying environment checksum."
        if [[ -e "$md5_file" ]]; then
            md5sum -c "$md5_file" >/dev/null 2>&1
            if [[ $? == 0 ]]; then
                install_needed=n
            fi
        fi

        if [[ $install_needed == y ]]; then
            read -p "Environment checksum did not match. Okay to delete and recreate $GGGPATH/install/.condaenv/? [yn]" answer
            case $answer in
                [yY])
                    echo "Removing $GGGPATH/install/.condaenv/"
                    rm -rf $GGGPATH/install/.condaenv/
                    ;;
                *)
                    echo "Not deleting the environment."
                    exit 2
            esac
        fi
fi

# I tried updating, but it's not reliable - it missed a toml dependency for some weird reason.
if [[ $install_needed == "y" ]]; then
        echo "Creating python environment from environment.yml"
        $CONDACMD $CREATECMD -p "$GGGPATH/install/.condaenv" -f environment.yml
        md5sum environment.yml > "$md5_file"
fi
