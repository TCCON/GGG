#!/bin/bash
always_yes=false

for arg in $@; do
    case $arg in 
        -y|--yes)
            always_yes=true
            ;;
    esac
done

srcdir="$GGGPATH/src/tccon_netcdf"
srcremote="https://bitbucket.org/rocheseb/tccon_netcdf"
#srcremote="https://jlaughner@bitbucket.org/jlaughner/tccon_netcdf.git"
#last_commit_for_ggg="c87fcf7"


printf "\n======== GETTING TCCON_NETCDF ========\n\n"

# In some cases, an active environment does not persist when a subshell
# is created by executing a script so we double check that here.
if [ $(basename $CONDA_PREFIX) != ggg-tccon-default ] ; then
    echo "The ggg-tccon-default python environment is not active. CONDA_PREFIX=$(basename ${CONDA_PREFIX}). Cannot proceed." 
    exit 1
fi

was_cloned=false
if ! [ -d $srcdir ] ; then
	git clone $srcremote $srcdir
	was_cloned=true
	cd $srcdir
else
	echo "$srcdir already exists"
fi

cd $srcdir
last_commit_for_ggg="a2ac07a"
git checkout master >/dev/null 2>/dev/null
git rev-parse --verify "ggg" &> /dev/null
if [ $? == 0 ] ; then
        if $always_yes; then
                answer=y
        else
                read -p "The ggg branch already exists in GGGPATH/src/tccon_netcdf; it will be overwritten and any new commits will be lost ; continue? [yn]: " answer
        fi

        case $answer in
                [yY])
                        echo "Proceeding ..."
                        ;;
                *)
                        echo "Aborting in clone_netcdf_writer.sh; if you did any work in the ggg branch of tccon_netcdf, move it elsewhere and rerun master.sh"
                        exit 1
                        ;;

        esac
fi
git branch -d ggg 2>/dev/null
echo "Creating ggg branch"
git checkout -b ggg
git reset --hard a803c94
git pull $srcremote master
git reset --hard $last_commit_for_ggg
git_head=$(git show --oneline -s | awk '{print $1}')
echo "The ggg branch HEAD is now at $git_head"

printf "\n======== INSTALLING TCCON_NETCDF ========\n\n"
./install_tccon_netcdf.sh || exit $?
printf "\n================= DONE =================\n\n"

