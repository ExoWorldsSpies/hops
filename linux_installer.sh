#!/bin/bash

mydir=$(pwd)

python3 "$mydir/setup.py" install

cp "$mydir/setup2_linux.py" "$HOME/setup2_linux.py"
python3 "$HOME/setup2_linux.py"
rm "$HOME/setup2_linux.py"
