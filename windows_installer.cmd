SET mydir=%cd%

python "%mydir%\setup.py" install

copy "%mydir%\setup2.py" "%HOMEDRIVE%%HOMEPATH%\setup2.py"
python "%HOMEDRIVE%%HOMEPATH%\setup2.py"
del "%HOMEDRIVE%%HOMEPATH%\setup2.py"