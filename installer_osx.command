#!/bin/bash

mydir="$(dirname "$BASH_SOURCE")"

python "$mydir/pre_setup.py"
python -m pip install -r "$mydir/requirements.txt"
python -m pip install "$mydir"
python "$mydir/post_setup.py"
python "$mydir/setup_close.py"
