#!/bin/bash

mydir="$(dirname "$BASH_SOURCE")"

cd "$mydir"

python "$mydir/pre_setup.py"
python -m pip install -r "$mydir/requirements.txt"
pip install -e .
python "$mydir/post_setup.py"
python "$mydir/setup_close.py"