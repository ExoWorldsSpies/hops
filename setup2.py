from __future__ import absolute_import
import platform
import os

from hops.pylightcurve3._1databases import Databases

databases = Databases(update=True)

from hops import __get_abspath__


name = 'hops'

app_dir = __get_abspath__()
system = platform.system()

executable = {'Darwin': 'command', 'Linux': 'sh', 'Windows': 'cmd'}

# create shortcut
try:
    shortcut = os.path.join(os.path.expanduser('~'), 'Desktop2', name + '.' + executable[system])
    w = open(shortcut, 'w')
    w.write('python \"' + app_dir + '\"')
    w.close()
except IOError:
    shortcut = os.path.join(os.path.expanduser('~'), name + '.' + executable[system])
    w = open(shortcut, 'w')
    w.write('python \"' + app_dir + '\"')
    w.close()
    print('\n\n************\nCould not find your Desktop. The shortcut is saved in your home folder:\n\n{0}\n\n'
          'You can freely move this file to your preferred location.\n************\n\n'.format(shortcut))
    x = input('Press enter to complete installation.\n')

if system == 'Darwin':
    os.system('chmod 755 ' + shortcut)
elif system == 'Linux':
    os.system('chmod +x ' + shortcut)
