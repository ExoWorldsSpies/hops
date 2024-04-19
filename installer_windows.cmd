SET mydir=%cd%

python "%mydir%\pre_setup.py"

SET /p pydir=<pydir.txt

call "%pydir%"

python -m pip install -r "%mydir%\requirements.txt"
python -m pip install "$mydir"
python "%mydir%\post_setup.py"
python "%mydir%\setup_close.py"
