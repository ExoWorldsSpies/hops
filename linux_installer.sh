#!/bin/bash

mydir=$(pwd)

python3 "$mydir/setup.py" install

cp "$mydir/setup2.py" "$HOME/setup2.py"
python3 "$HOME/setup2_linux.py"
rm "$HOME/setup2.py"
