
__all__ = ['HOPS']

import os
import glob
import datetime
import webbrowser
import hops.pylightcurve41 as plc

from .application_log import HOPSLog
from .application_windows import MainWindow
from .application_1_data_target import DataTargetWindow
from .application_2_reduction import ReductionWindow
from .application_3_inspectiont import InspectiontWindow
from .application_4_alignment import AlignmentWindow
from .application_5_photometry import PhotometryWindow, PhotometryProgressWindow
from .application_6_fitting import FittingWindow
from .application_extra_1_observing_planner import ObservingPlannerWindow


class HOPS(MainWindow):

    def __init__(self):

        MainWindow.__init__(self, HOPSLog(), name='HOlomon Photometric Software', position=1)

        # logo

        logo = self.FigureWindow(figsize=(1.8, 1.8))

        ax = logo.figure.add_subplot(111)
        ax.imshow(self.log.logo_jpg)
        ax.axis('off')
        logo.figure.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

        self.my_profile_window = None
        self.my_profile_core_headers = None
        self.my_profile_labels = None
        self.my_profile_entries = None

        # directory

        self.data_target_button = self.Button(text='** SELECT DATA & TARGET **', command=self.select_data_target)
        self.data_target_complete = self.Label()

        # reduction

        self.reduction_button = self.Button(text='** RUN REDUCTION **', command=self.run_reduction)
        self.reduction_complete = self.Label()

        # for backwards compatibility
        if self.log.get_param('reduction_complete'):
            if not self.log.is_version_new(self.log.get_param('reduction_version'), '3.0.0'):
                self.log.set_param('reduction_complete', False)

        # inspection

        self.inspection_button = self.Button(text='INSPECT FRAMES', command=self.inspect_frames)
        self.inspection_complete = self.Label()

        # alignment

        self.alignment_button = self.Button(text='** RUN ALIGNMENT **', command=self.run_alignment)
        self.alignment_complete = self.Label()

        # for backwards compatibility
        if self.log.get_param('alignment_complete'):
            if not self.log.is_version_new(self.log.get_param('reduction_version'), '3.0.0'):
                self.log.set_param('alignment_complete', False)

        # photometry

        self.photometry_button = self.Button(text='** PHOTOMETRY **', command=self.photometry)
        self.photometry_complete = self.Label()

        # for backwards compatibility
        if self.log.get_param('photometry_complete'):
            if not self.log.is_version_new(self.log.get_param('photometry_version'), '3.0.0'):
                self.log.set_param('photometry_complete', False)

        # fitting

        self.fitting_button = self.Button(text='EXOPLANET FITTING', command=self.fitting)
        self.fitting_complete = self.Label()

        if self.log.get_param('fitting_complete'):
            if not self.log.is_version_new(self.log.get_param('fitting_version'), '3.0.0'):
                self.log.set_param('fitting_complete', False)

        #

        self.setup_window([
            [[logo, 0, 2, 6]],
            [[self.Label(text='{0}\nv{1}'.format(self.log.software_name, self.log.version)), 2, 1, 1, 'title']],
            [[self.Label(text='Copyright (c) 2017-{1} Angelos Tsiaras, atsiaras@star.ucl.ac.uk'.format(
              self.log.version, datetime.date.today().year)), 2]],
            [[self.Button(text=self.log.updates, command=self.open_updates), 2, 1, 2]],
            [],
            [[self.Button(text='MY PROFILE', command=self.open_my_profile), 2]],
            [],
            [[self.Label(text='Analyse your data step by step'), 0, 3, 1, 'title']],
            [],
            [[self.Label(text='     1. '), 0], [self.data_target_button, 1], [self.data_target_complete, 2]],
            [],
            [[self.Label(text='     2. '), 0], [self.reduction_button, 1], [self.reduction_complete, 2]],
            [],
            [[self.Label(text='     3. '), 0], [self.inspection_button, 1], [self.inspection_complete, 2]],
            [],
            [[self.Label(text='     4. '), 0], [self.alignment_button, 1], [self.alignment_complete, 2]],
            [],
            [[self.Label(text='     5. '), 0], [self.photometry_button, 1], [self.photometry_complete, 2]],
            [],
            [[self.Label(text='     6. '), 0], [self.fitting_button, 1], [self.fitting_complete, 2]],
            [],
            [[self.Label(text='** mandatory step **'), 1], [self.Button(text='EXIT', command=self.close), 2]],
            [],
            [[self.Label(text='Extra tools:'), 1]],
            [[self.Button(text='OBSERVING PLANNER', command=self.open_observing_planner), 1]],
            []
        ])

        self.update_window()

    def open_updates(self):
            webbrowser.open("https://www.exoworldsspies.com/en/software/", new=1)

    def open_my_profile(self):

        self.my_profile_window = MainWindow(self.log, name='HOPS - My Profile', position=1)

        # my profile window

        self.my_profile_core_headers = self.log.main_log_profile.keys()

        self.my_profile_labels = {}
        self.my_profile_entries = {}
        for row, header in enumerate(self.my_profile_core_headers):
            self.my_profile_labels[header] = self.my_profile_window.Label(text=header)
            try:
                self.my_profile_entries[header] = self.my_profile_window.Entry(value=self.log.local_log_profile[header])
            except:
                self.my_profile_entries[header] = self.my_profile_window.Entry(value=self.log.main_log_profile[header])

        elements_per_column = 17
        structure = [
            [[self.my_profile_window.Label(text='Header Keywords for useful information,\nif included in the FITS header:'), 0, 2],
             [self.my_profile_window.Label(text='Default values for useful information,\nif not included in the FITS header:'), 2, 2],
             [self.my_profile_window.Label(text='\nDefault file names:'), 4, 2]
             ],
            [],
        ]

        for ff in range(elements_per_column):
            structure.append([])

        for jj, header in enumerate(list(self.my_profile_core_headers)):

            column = int(jj/elements_per_column)
            row = jj - column * elements_per_column + 2

            if '___' in header:
                if column == 1:
                    structure[row].append([self.my_profile_window.Label(text='   Not Applicable   '), column * 2, 2])
                else:
                    pass
            else:
                structure[row].append([self.my_profile_labels[header], column * 2])
                structure[row].append([self.my_profile_entries[header], column * 2 + 1])

        structure[8].append([self.my_profile_window.Label(text='Information used only by the scheduler:'), 4, 2])

        structure.append([])

        structure.append([[self.my_profile_window.Button(text='SAVE CHANGES & CLOSE WINDOW',
                                                         command=self.update_local_log_profile), 0, 6]])
        structure.append([])

        self.my_profile_window.setup_window(structure)

        self.my_profile_window.run()

    def update_local_log_profile(self):
        test = {}
        for header in self.my_profile_core_headers:
            test[header] = self.my_profile_entries[header].get()
        self.log.save_local_log_profile(test)

        self.my_profile_window.close()

    def select_data_target(self):

        self.log.set_param('proceed', False)

        self.disable()

        data_target_window = DataTargetWindow(self.log)
        data_target_window.run(f_after=data_target_window.check_dir)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            self.run_reduction()

    def run_reduction(self):

        self.log.set_param('proceed', False)

        self.disable()

        reduction_window = ReductionWindow(self.log)
        reduction_window.run(f_before=reduction_window.progress_figure.adjust_size,
                             f_after=reduction_window.run_reduction)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            self.inspect_frames()

    def inspect_frames(self):

        self.log.set_param('proceed', False)

        self.disable()

        inspection_window = InspectiontWindow(self.log)
        inspection_window.run(f_before=inspection_window.adjust_size)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            self.run_alignment()

    def run_alignment(self):

        self.log.set_param('proceed', False)

        self.disable()

        alignment_window = AlignmentWindow(self.log)
        alignment_window.run(f_before=alignment_window.progress_figure.adjust_size,
                             f_after=alignment_window.run_alignment)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            self.photometry()

    def photometry(self):

        self.log.set_param('proceed', False)

        self.disable()

        photometry_window = PhotometryWindow(self.log)
        photometry_window.run(f_before=photometry_window.fits_figure.adjust_size)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            if self.log.get_param('proceed') == 'run_photometry':
                self.run_photometry()
            else:
                self.fitting()

    def run_photometry(self):

        self.log.set_param('proceed', False)

        self.disable()

        photometry_progress_window = PhotometryProgressWindow(self.log)
        photometry_progress_window.run(
            f_before=photometry_progress_window.progress_figure.adjust_size,
            f_after=photometry_progress_window.run_photometry)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            if self.log.get_param('proceed') == 'return_to_photometry':
                self.photometry()
            else:
                self.fitting()

    def fitting(self):

        self.log.set_param('proceed', False)

        self.disable()

        fitting_window = FittingWindow(self.log)
        fitting_window.run(f_after=fitting_window.preview_figure.adjust_size)

        self.activate()
        self.update_window()

        if self.log.get_param('proceed'):
            if self.log.get_param('proceed') == 'return_to_photometry':
                self.photometry()

    def open_observing_planner(self):

        self.disable()

        observing_planner_window = ObservingPlannerWindow(self.log)
        observing_planner_window.run()

        self.activate()

    def update_window(self):

        self.data_target_button.activate()

        if self.log.get_param('data_target_complete'):

            self.data_target_complete.set('Data:     {0}\nTarget:    {1} '.format(
                self.log.get_param('directory_short'), self.log.get_param('target_name')))

            self.reduction_button.activate()
            trash = 0

            if not os.path.isfile(self.log.all_frames):
                self.log.set_param('reduction_complete', False)
            else:
                all_frames = plc.open_dict(self.log.all_frames)
                for i in all_frames:
                    if not os.path.isfile(os.path.join(self.log.reduction_directory, i)):
                        self.log.set_param('reduction_complete', False)
                        break
                    if all_frames[i][self.log.skip_key]:
                        trash += 1

            if self.log.get_param('reduction_complete'):
                self.reduction_complete.set('Completed under v{0}'.format(self.log.get_param('reduction_version')))
                self.inspection_button.activate()
                self.inspection_complete.set('Files discarded: {0}'.format(trash))
                self.alignment_button.activate()

                all_frames = plc.open_dict(self.log.all_frames)
                for i in all_frames:
                    if not all_frames[i][self.log.skip_key] and not all_frames[i][self.log.align_x0_key]:
                        self.log.set_param('alignment_complete', False)
                plc.save_dict(all_frames, self.log.all_frames)

                if not os.path.isfile(self.log.all_stars):
                    self.log.set_param('alignment_complete', False)

                if self.log.get_param('alignment_complete'):

                    self.alignment_complete.set('Completed under v{0}'.format(self.log.get_param('alignment_version')))
                    self.photometry_button.activate()

                    if len(glob.glob('{0}*'.format(self.log.photometry_directory_base))) == 0:
                        self.log.set_param('photometry_complete', False)

                    if self.log.get_param('photometry_complete'):
                        self.photometry_complete.set('Completed under v{0}'.format(self.log.get_param('photometry_version')))
                        self.fitting_button.activate()

                        if self.log.get_param('fitting_complete'):
                            self.fitting_complete.set('Completed under v{0}'.format(self.log.get_param('fitting_version')))

                    else:
                        self.photometry_complete.set('You need to complete this step to proceed')
                        self.fitting_button.disable()
                        self.fitting_complete.set('')


                else:
                    self.alignment_complete.set('You need to complete this step to proceed')
                    self.photometry_button.disable()
                    self.photometry_complete.set('')
                    self.fitting_button.disable()
                    self.fitting_complete.set('')

            else:
                self.reduction_complete.set('You need to complete this step to proceed')
                self.inspection_button.disable()
                self.inspection_complete.set('')
                self.alignment_button.disable()
                self.alignment_complete.set('')
                self.photometry_button.disable()
                self.photometry_complete.set('')
                self.fitting_button.disable()
                self.fitting_complete.set('')

        else:
            self.data_target_complete.set('You need to complete this step to proceed.')
            self.reduction_button.disable()
            self.reduction_complete.set('')
            self.inspection_button.disable()
            self.inspection_complete.set('')
            self.alignment_button.disable()
            self.alignment_complete.set('')
            self.photometry_button.disable()
            self.photometry_complete.set('')
            self.fitting_button.disable()
            self.fitting_complete.set('')
