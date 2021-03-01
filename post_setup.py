
import os
import glob
import platform
import sys
sys.path = sys.path[1:]

from hops import __get_abspath__ as hops__get_abspath__

name = 'hops'

here = os.path.abspath(os.path.dirname(__file__))
shortcut = open(os.path.join(here, 'shortcut.txt')).read()
shortcut = shortcut.replace('{{x}}', hops__get_abspath__())

# create shortcut
print('\n\n************\n')

executable = {'Darwin': 'command', 'Linux': 'sh', 'Windows': 'cmd'}[platform.system()]

try:
    shortcut_path = os.path.join(os.path.expanduser('~'), 'Desktop', name + '.' + executable)
except IOError:
    try:
        shortcut_path = os.path.join(glob.glob(os.path.join(os.path.expanduser('~'), '*', 'Desktop'))[0], name + '.' + executable)
    except IndexError:
        shortcut_path = os.path.join(os.path.expanduser('~'), name + '.' + executable)

w = open(shortcut_path, 'w')
w.write(shortcut)
w.close()

if platform.system() == 'Darwin':
    os.system('chmod 755 ' + shortcut_path)
elif platform.system() == 'Linux':
    os.system('chmod +x ' + shortcut_path)


print('HOPS successfully installed.')
print('The shortcut has been saved in your Desktop:\n\n{0}\n\n'
      'You can freely move this file to your preferred location.'.format(shortcut_path))
