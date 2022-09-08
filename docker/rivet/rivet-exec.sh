#! /usr/bin/env bash
# Wrap the bash setup so it'll still work as a Singularity container

. /etc/profile
. /etc/bash.bashrc
for i in /etc/profile.d/*; do source $i; done

# export PYTHONPATH=/usr/local/lib/python3.8/site-packages:/usr/local/lib/python3.8/dist-packages:$PYTHONPATH
# export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
# echo "::$PYTHONPATH::"
# ls /usr/local/lib/python3.8/site-packages
# python -c "import rivet"

if [[ -n "$@" ]]; then
    # python -c "import sys"
    # python -c "import commands"
    # python -c "import rivet"
    eval "$@"
else
    bash --rcfile /etc/bash.bashrc
fi
