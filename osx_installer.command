#!/bin/bash

mydir="$(dirname "$BASH_SOURCE")"

python "$mydir/setup.py" install

cp "$mydir/setup2.py" "$HOME/setup2.py"
python "$HOME/setup2.py"
rm "$HOME/setup2.py"