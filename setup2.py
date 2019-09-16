from __future__ import absolute_import
import platform
import os
from hops import __get_abspath__

name = 'hops'

app_dir = __get_abspath__()
system = platform.system()

excecutable = {'Darwin': 'command', 'Linux': 'sh', 'Windows': 'cmd'}

# create shortcut
try:
    shortcut = os.path.join(os.path.expanduser('~'), 'Desktop', name + '.' + excecutable[system])
    w = open(shortcut, 'w')
    w.write('python \"' + app_dir + '\"')
    w.close()
except IOError:
    shortcut = os.path.join(os.path.expanduser('~'), name + '.' + excecutable[system])
    w = open(shortcut, 'w')
    w.write('python \"' + app_dir + '\"')
    w.close()

if system == 'Darwin':
    os.system('chmod 755 ' + shortcut)
elif system == 'Linux':
    os.system('chmod +x ' + shortcut)
