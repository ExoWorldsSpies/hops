#!/bin/bash

mydir="$(dirname "$BASH_SOURCE")"

python "$mydir/pre_setup.py"
python -m pip install -r requirements.txt
python "$mydir/setup.py" install
python "$mydir/post_setup.py"
python "$mydir/setup_close.py"
