from __future__ import absolute_import
import platform
import os
import glob

from hops.pylightcurve3._1databases import Databases

databases = Databases(update=True)

from hops import __get_abspath__


name = 'hops'

app_dir = __get_abspath__()
system = platform.system()

executable = {'Darwin': 'command', 'Linux': 'sh', 'Windows': 'cmd'}

# create shortcut
print('\n\n************\n')
try:
    shortcut = os.path.join(os.path.expanduser('~'), 'Desktop', name + '.' + executable[system])
    w = open(shortcut, 'w')
    w.write('python \"' + app_dir + '\"')
    w.close()
    print('HOPS successfully installed.')
    print('The shortcut has been saved in your Desktop:\n\n{0}\n\n'
          'You can freely move this file to your preferred location.'.format(shortcut))
except IOError:
    try:
        shortcut = os.path.join(glob.glob(os.path.join(os.path.expanduser('~'), '*', 'Desktop'))[0], name + '.' + executable[system])
        w = open(shortcut, 'w')
        w.write('python \"' + app_dir + '\"')
        w.close()
        print('HOPS successfully installed.')
        print('The shortcut has been saved in your Desktop:\n\n{0}\n\n'
              'You can freely move this file to your preferred location.'.format(shortcut))
    except IndexError:
        shortcut = os.path.join(os.path.expanduser('~'), name + '.' + executable[system])
        w = open(shortcut, 'w')
        w.write('python \"' + app_dir + '\"')
        w.close()
        print('HOPS successfully installed.')
        print('You Desktop could not be located. The shortcut has been saved in your home folder:\n\n{0}\n\n'
              'You can freely move this file to your preferred location.'.format(shortcut))


if system == 'Darwin':
    os.system('chmod 755 ' + shortcut)
elif system == 'Linux':
    os.system('chmod +x ' + shortcut)


print('\n************\n\n')
x = input('Press enter to exit.\n')