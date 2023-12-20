#!/bin/bash

mydir=$(pwd)

python3 "$mydir/pre_setup.py"
python3 -m pip install -r requirements.txt
python3 "$mydir/setup.py" install
python3 "$mydir/post_setup.py"
python3 "$mydir/setup_close.py"
