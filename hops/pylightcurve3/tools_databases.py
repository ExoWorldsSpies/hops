from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._0errors import *

import os
import sys
import glob
import time
import pickle
import shutil
import datetime
import ssl
import urllib
from urllib.request import urlretrieve

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


class Database:

    def __init__(self, database_name, package_name, file=__file__, vital=False, date_to_update='daily', expire_date=None,
                 frequencey=None, force_update=False, ask_size=None, database_file_name='_0database.pickle',
                 package_path=os.path.expanduser('~'), directory_name='database',
                 last_update_file_name='database_last_update.txt'):

        info_file_path = os.path.join(os.path.abspath(os.path.dirname(file)), database_file_name)
        package_path = os.path.join(package_path, '.{0}'.format(package_name))
        if not os.path.isdir(package_path):
            os.mkdir(package_path)
        directory_path = os.path.join(package_path, '{0}_{1}'.format(database_name, directory_name))
        last_update_file_path = os.path.join(package_path, '{0}_{1}'.format(database_name, last_update_file_name))

        if os.path.isdir(directory_path):
            if force_update or len(glob.glob(os.path.join(directory_path, '*'))) == 0:
                shutil.rmtree(directory_path)
                os.mkdir(directory_path)
                update = True
            else:
                if not os.path.isfile(last_update_file_path):
                    update = True

                else:

                    last_update = open(last_update_file_path).readlines()[0]
                    last_update_date = datetime.date(int('20{0}'.format(last_update[:2])), int(last_update[2:4]),
                                                     int(last_update[4:6]))
                    must_update_date = datetime.date(int('20{0}'.format(date_to_update[:2])), int(date_to_update[2:4]),
                                                     int(date_to_update[4:6]))

                    if frequencey:
                        must_update_date_2 = last_update_date + datetime.timedelta(days=frequencey)
                    else:
                        must_update_date_2 = datetime.datetime.now().date() + datetime.timedelta(days=365)

                    if expire_date:
                        must_update_date_3 = datetime.date(int('20{0}'.format(expire_date[:2])), int(expire_date[2:4]),
                                                           int(expire_date[4:6]))
                    else:
                        must_update_date_3 = must_update_date_2

                    if (last_update_date < must_update_date or
                            datetime.datetime.now().date() > must_update_date_2 or
                            datetime.datetime.now().date() > must_update_date_3):
                        update = True
                    else:
                        update = False
        else:
            os.mkdir(directory_path)
            update = True

        if update and ask_size:
            if input('Downloading {0} database (up to {1})... proceed with download now? (y/n): '.format(
                    database_name, ask_size)) == 'y':
                update = True
            else:
                update = False

        if update:
            # noinspection PyBroadException
            try:
                print('\nDownloading {0} database...'.format(database_name))

                dbx_files = pickle.load(open(info_file_path, 'rb'))
                dbx_files = dbx_files['{0}_{1}'.format(database_name, directory_name)]

                if database_name == 'clablimb':
                    if len(glob.glob(os.path.join(directory_path, '*.pickle'))) > 0:
                        if os.path.isfile(glob.glob(os.path.join(directory_path, '*.pickle'))[0]):
                            try:
                                xx = pickle.load(open(glob.glob(os.path.join(directory_path, '*.pickle'))[0], 'rb'))
                                del xx
                            except:
                                os.remove(glob.glob(os.path.join(directory_path, '*.pickle'))[0])

                for current_file in glob.glob(os.path.join(directory_path, '*')):
                    if os.path.split(current_file)[1] not in dbx_files:
                        os.remove(current_file)

                for dbx_file in dbx_files:
                    if not os.path.isfile(os.path.join(package_path, dbx_files[dbx_file]['local_path'])):
                        print(dbx_file)
                        try:
                            urlretrieve(dbx_files[dbx_file]['link'], os.path.join(package_path,
                                                                                  dbx_files[dbx_file]['local_path']))
                        except:
                            with urllib.request.urlopen(dbx_files[dbx_file]['link'], context=ctx) as u, \
                                        open(os.path.join(package_path, dbx_files[dbx_file]['local_path']), 'wb') as f:
                                f.write(u.read())

                for update_file in glob.glob(os.path.join(directory_path, 'update_*')):
                    shutil.move(update_file, update_file.replace('update_', ''))

                if database_name == 'clablimb':
                    xx = pickle.load(open(glob.glob(os.path.join(directory_path, '*.pickle'))[0], 'rb'))
                    for i in xx:
                        w = open(os.path.join(directory_path, i), 'w')
                        w.write(xx[i])
                        w.close()

                w = open(last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

            except Exception as inst:
                print('\nDownloading {0} database failed. A download will be attempted next time.'.format(
                    database_name))
                print('Error:', sys.exc_info()[0])
                print(inst.args)
                pass

        if (not os.path.isdir(directory_path) or
                len(glob.glob(os.path.join(directory_path, '*'))) == 0):
            if vital:
                raise PyLCLibraryError('{0} database not available.'.format(database_name))
            else:
                print('\n{0} features cannot be used.'.format(database_name))
                self.path = False
        else:
            self.path = directory_path
