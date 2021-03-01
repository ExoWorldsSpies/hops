
import sys
import os
import platform
import time

system = platform.system()
python_dir = os.path.split(sys.executable)[0]

here = os.path.abspath(os.path.dirname(__file__))

w = open(os.path.join(here, 'shortcut.txt'), 'w')
if platform.system() == 'Windows':
    if 'anaconda' in sys.executable.lower():
        if os.path.isfile(os.path.join(python_dir, 'Scripts', 'activate.bat')):
            w.write('call "' + os.path.join(python_dir, 'Scripts', 'activate.bat') + '"\npython {{x}}')
            os.system('echo {0}>pydir.txt'.format(os.path.join(python_dir, 'Scripts', 'activate.bat')))
        else:
            w.write('python {{x}}')
            os.system('echo {0}>pydir.txt'.format(os.path.join(here, 'activate.bat')))
    else:
        w.write('python {{x}}')
        os.system('echo {0}>pydir.txt'.format(os.path.join(here, 'activate.bat')))

elif platform.system() == 'Linux':
    w.write('python3 {{x}}')
else:
    w.write('python {{x}}')

w.close()

time.sleep(1)
