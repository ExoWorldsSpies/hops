SET mydir=%cd%

python "%mydir%\pre_setup.py"

SET /p pydir=<pydir.txt

call "%pydir%"

python -m pip install -r requirements.txt
python "%mydir%\setup.py" install
python "%mydir%\post_setup.py"
python "%mydir%\setup_close.py"
