#!/bin/bash

cd "$GGGPATH/install/"
source "bash_install_utils.sh"

# Necessary to use conda commands inside this script, even though this was already 
# done in run_pyinstall.sh. This will set a bunch of variables that we use to run the
# installation. If you can't find a variable in the code, check this file (created by
# check_python.sh)
source ".init_conda"

md5_file="${ENVDIR}/environment.md5"

verify_environment() {
    env_dir="$1"
    install_needed=y
    if [[ -e "$env_dir" ]]; then
            echo "The $env_dir python environment already exists; verifying environment checksum."
            if [[ -e "$md5_file" ]]; then
                md5sum -c "$md5_file" >/dev/null 2>&1
                if [[ $? == 0 ]]; then
                    install_needed=n
                fi
            fi

            if [[ $install_needed == y ]]; then
                # I tried updating the conda environment in place, but it's not reliable.
                # It missed a toml dependency for some weird reason.
                read -p "Environment checksum did not match. Okay to delete and recreate ${env_dir}? [yn]" answer
                case $answer in
                    [yY])
                        echo "Removing ${env_dir}"
                        rm -rf "$env_dir"
                        ;;
                    *)
                        echo "Not deleting the environment."
                        return 2
                esac
            fi
    fi

    if [[ "$install_needed" == "n" ]]; then
        # I picked a high number unlikely to be a normal exit
        # code but which still fits in an 8 bit signed or unsiged
        # integer, in case different systems use different conventions
        # for return codes...
        return 100
    else
        return 0
    fi
}



verify_environment "$ENVDIR"
rcode=$?
if [[ $rcode == 100 ]]; then
    # Indicates installation is not required
    echo "Environment is up to date, not installing."
    exit 0
elif [[ $rcode -gt 0 ]]; then
    exit $rcode
fi

if do_use_pip_only; then
    echo "Creating python environment from requirements.txt at ${ENVDIR}"
    "$PY_FOR_CREATE" -m venv "$ENVDIR"
    source "$ENVDIR/bin/activate"
    "$PIPCMD" install --requirement requirements.txt
    md5sum requirements.txt > "$md5_file"
else
    echo "Creating python environment from environment.yml at ${ENVDIR}"
    $CONDACMD $CREATECMD -p "$ENVDIR" -f environment.yml
    md5sum environment.yml > "$md5_file"
fi
