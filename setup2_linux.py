from __future__ import absolute_import
import os
from hops import __get_abspath__

name = 'hops'

app_dir = __get_abspath__()

# create shortcut
try:
    shortcut = os.path.join(os.path.expanduser('~'), 'Desktop', name + '.sh')
    w = open(shortcut, 'w')
    w.write('python3 \"' + app_dir + '\"')
    w.close()
except IOError:
    shortcut = os.path.join(os.path.expanduser('~'), name + '.sh')
    w = open(shortcut, 'w')
    w.write('python3 \"' + app_dir + '\"')
    w.close()

os.system('chmod +x ' + shortcut)
