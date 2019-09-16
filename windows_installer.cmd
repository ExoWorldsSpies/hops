SET mydir=%cd%

python "%mydir%\setup.py" install

copy "%mydir%\setup2.py" "%HOMEPATH%\setup2.py"
python "%HOMEPATH%\setup2.py"
del "%HOMEPATH%\setup2.py"