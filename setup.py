
import os
import glob

from setuptools import setup

name = 'hops'
description = 'HOlomon Photometry Software - A software to analyse data from small ground-based telescopes'
url = 'https://https://github.com/ExoWorldsSpies/hops'
install_requires = ['pyaml', 'requests', 'matplotlib>=3.3.3', 'numpy>=1.20.0', 'emcee', 'seaborn',
                    'astropy', 'scipy', 'astroquery>=0.4', 'pillow', 'quantities']

os.chdir(os.path.abspath(os.path.dirname(__file__)))

subdirs_to_include = []
for x in os.walk(name):
    if os.path.isdir(x[0]):
        if x[0] != name:
            subdirs_to_include.append(x[0])

files_to_include = []
for x in glob.glob(os.path.join(name, '*')):
    if x[-2:] != 'py':
        files_to_include.append(os.path.join(name, os.path.split(x)[1]))

files_to_include.append('README.md')
files_to_include.append('LICENSE')

w = open('MANIFEST.in', 'w')
for i in subdirs_to_include:
    w.write('include ' + os.path.join(i, '*') + ' \n')

for i in files_to_include:
    w.write('include ' + i + ' \n')

w.close()

version = ' '
for i in open(os.path.join(name, '__init__.py')):
    if len(i.split('__version__')) > 1:
        version = i.split()[-1][1:-1]

print('HOPS ', version)

setup(
    name=name,
    version=version,
    description=description,
    url=url,
    author='Angelos Tsiaras',
    author_email='aggelostsiaras@gmail.com',
    license='MIT',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X'
                 'Programming Language :: Python :: 3.7',
                 ],
    packages=[name],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)
