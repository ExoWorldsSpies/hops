
__all__ = ['MainWindow', 'AddOnWindow']

import os
import time
import numpy as np
import datetime
import matplotlib
import webbrowser
import warnings

from tkinter import Tk, TclError
from tkinter import Label, Button, Entry, Checkbutton, Scrollbar, Listbox, PhotoImage, Radiobutton, Scale, Frame, Canvas
from tkinter import StringVar, BooleanVar, DoubleVar, IntVar
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, TOP, BOTTOM, BOTH, Y, HORIZONTAL, VERTICAL, E, W, N, S, NW, TRUE, FALSE

from tkinter.ttk import Combobox, Style, Progressbar
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showinfo, askyesno, askyesnocancel

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    NavigationToolbar2Tk = NavigationToolbar2TkAgg

from matplotlib.cm import Greys, Greys_r
from astropy.io import fits as pf

import hops.pylightcurve41 as plc
from hops.hops_tools.fits import get_fits_data
from hops.application_log import HOPSLog


matplotlib.use('TkAgg')


def openweb():
    webbrowser.open("https://www.exoworldsspies.com/en/software", new=1)


class HOPSWindow:

    def __init__(self, log=HOPSLog, name='HOPS - window', sizex=None, sizey=None, position=5):

        self.log = log

        self.root = Tk()
        self.hide()
        self.main_frame = self.root
        self.root.protocol('WM_DELETE_WINDOW', self.close)

        self.name = name
        self.position = position
        self.root.wm_title(name)
        if sizex and sizey:
            self.root.geometry('{0}x{1}'.format(int(self.root.winfo_screenwidth() / sizex),
                                                int(self.root.winfo_screenheight() / sizey)))

        self.mainloop_on = False
        self.exit = False
        self.widgets = []
        self.jobs = []
        self.jobs_completed = 0

        self.DISABLED = DISABLED
        self.NORMAL = NORMAL
        self.END = END
        self.RIGHT = RIGHT
        self.LEFT = LEFT
        self.TOP = TOP
        self.BOTTOM = BOTTOM
        self.BOTH = BOTH
        self.Y = Y
        self.HORIZONTAL = HORIZONTAL
        self.VERTICAL = VERTICAL

    def no_action(self):
        pass

    def update(self):

        self.root.update()

    def after(self, function=None, time=10):

        if function:

            def nested_after(function1, function2):

                def xx():
                    function1()
                    self.root.after(time, function2)

                return xx

            def internal_command():

                if isinstance(function, list):
                    functions = [nested_after(function[-2], function[-1])]
                    for i in range(len(function) - 2):
                        functions.append(nested_after(function[-i-3], functions[-1]))
                    self.root.after(time, functions[-1])

                else:
                    function()

                self.jobs_completed += 1

            xx = self.root.after(time, internal_command)
            self.jobs.append(xx)
            self.update_idletasks()

    def update_idletasks(self):

        self.root.update_idletasks()

    def reposition(self):

        self.root.update_idletasks()

        if self.position == 1:
            x = 0
            y = 0

        elif self.position == 2:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = 0

        elif self.position == 3:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = 0

        elif self.position == 4:
            x = 0
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 5:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 6:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 7:
            x = 0
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        elif self.position == 8:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        elif self.position == 9:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        elif self.position == 10:
            x = self.root.winfo_screenwidth()/2 - self.root.winfo_reqwidth()
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 11:
            x = self.root.winfo_screenwidth()/2
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        else:
            x = 0
            y = 0

        self.root.geometry('+%d+%d' % (int(x), int(y)))

        self.root.update_idletasks()

    def show(self):

        self.reposition()

        self.root.wm_attributes("-topmost", 1)
        self.root.after_idle(self.root.attributes, '-topmost', 0)

        self.root.deiconify()
        self.update_idletasks()

    def hide(self):

        self.root.withdraw()

    def run(self, f_after=None, f_before=None):

        print('\nStarting window: ', self.name)

        if f_before:
            f_before()

        self.show()

        self.mainloop_on = True

        self.after(f_after)

        self.root.mainloop()

    def trigger_exit(self):
        self.exit = True

    def def_close(self):

        if self.mainloop_on:
            self.root.quit()

        for job in self.jobs:
            self.root.after_cancel(job)

        self.root.destroy()

    def close(self):

        self.def_close()

    # widgets

    def register(self, widget):

        self.widgets.append(widget)

    def disable(self):

        self.root.protocol('WM_DELETE_WINDOW', self.no_action)

        for widget in self.widgets:
            widget.disable()

    def activate(self):

        self.root.protocol('WM_DELETE_WINDOW', self.close)

        for widget in self.widgets:
            widget.activate()

    # help pop ups

    def askdirectory(self):

        return askdirectory()

    def askyesno(self, *args, **kwargs):

        return askyesno(*args, **kwargs)

    def askyesnocancel(self, *args, **kwargs):

        return askyesnocancel(*args, **kwargs)

    def showinfo(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            showinfo(*args, **kwargs)

    # tk variables initiated directly in - and linked with - the hops window

    def StringVar(self, value):
        return StringVar(self.root, value=value)

    def BooleanVar(self, value):
        return BooleanVar(self.root, value=value)

    def DoubleVar(self, value):
        return DoubleVar(self.root, value=value)

    def IntVar(self, value):
        return IntVar(self.root, value=value)

    # tk widgwts initiated directly in - and linked with - the hops window

    def Label(self, **kwargs):
        return HOPSLabel(self, **kwargs)

    def Entry(self, **kwargs):
        return HOPSEntry(self, **kwargs)

    def Button(self, **kwargs):
        return HOPSButton(self, **kwargs)

    def DropDown(self, **kwargs):
        return HOPSDropDown(self, **kwargs)

    def CheckButton(self, **kwargs):
        return HOPSCheckButton(self, **kwargs)

    def ListDisplay(self, **kwargs):
        return HOPSListDisplay(self,  **kwargs)

    def Progressbar(self, **kwargs):
        return HOPSProgressbar(self, **kwargs)

    def FigureWindow(self, **kwargs):
        return HOPSFigureWindow(self, **kwargs)

    def FitsWindow(self, **kwargs):
        return HOPSFitsWindow(self, **kwargs)

    def Frame(self, *args, **kwargs):
        return Frame(self.main_frame, *args, **kwargs)

    def Radiobutton(self, *args, **kwargs):
        return Radiobutton(self.main_frame, *args, **kwargs)

    def setup_window(self, objects, title_font='times', main_font='times', button_font='times', entries_wd=20, entries_bd=3, buttons_bd=5):

        ff = min(1, (self.root.winfo_screenheight() / self.root.winfo_fpixels('1i')) / 13)
        font_size = 15 * ff

        button_font = (button_font, int(1.1 * font_size), 'bold')
        main_font = (main_font, int(font_size))
        title_font = (title_font, int(1.5 * font_size), 'bold')

        for row in range(len(objects)):
            if len(objects[row]) == 0:
                label_empty = Label(self.main_frame, text='')
                label_empty.grid(row=row, column=100)
            else:
                for obj in objects[row]:

                    if obj[0].winfo_class() == 'Button':
                        obj[0].config(borderwidth=buttons_bd, font=button_font, padx=1, pady=1)
                    elif obj[0].winfo_class() == 'Entry':
                        obj[0].configure(width=entries_wd, bd=entries_bd, font=main_font)
                    elif obj[0].winfo_class() in ['Label', 'Radiobutton', 'Checkbutton']:
                        if len(obj) == 5:
                            if obj[4] == 'title':
                                obj[0].configure(font=title_font, padx=0, pady=0)
                            else:
                                obj[0].configure(font=main_font, padx=0, pady=0)
                        else:
                            obj[0].configure(font=main_font, padx=0, pady=0)

                    if len(obj) >= 4:
                        obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                    elif len(obj) == 3:
                        obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                    else:
                        obj[0].grid(row=row, column=obj[1])


class MainWindow(HOPSWindow):

    def __init__(self, log, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self, log, name, sizex, sizey, position)

        self.registered_addon_windows = []

    def register_addon_window(self, addon_window):
        self.registered_addon_windows.append(addon_window)

    def close(self):

        for addon_window in self.registered_addon_windows:
            addon_window.def_close()

        self.def_close()

    def set_close_button_function(self, function):
        self.root.protocol('WM_DELETE_WINDOW', function)


class AddOnWindow(HOPSWindow):

    def __init__(self, main_window, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self, main_window.log, name, sizex, sizey, position)

        main_window.register_addon_window(self)

    def close(self):
        self.hide()

# widgets


class HOPSWidget:

    def __init__(self, window, widget, name):

        self.window = window
        window.register(self)
        self.widget = widget
        self.name = name
        self.root = self.window.root
        self.log = self.window.log

    def winfo_class(self):
        return self.name

    def grid(self, *args, **kwargs):
        self.widget.grid(*args, **kwargs)

    def disable(self):
        if self.name not in ['FigureWindow', 'Progressbar']:
            self.widget['state'] = DISABLED

    def activate(self):
        if self.name in ['FigureWindow', 'Progressbar']:
            pass
        elif self.name == 'DropDown':
            self.widget['state'] = 'readonly'
        else:
            self.widget['state'] = NORMAL

    def configure(self,  *args, **kwargs):
        self.widget.configure(*args, **kwargs)

    def config(self,  *args, **kwargs):
        self.widget.configure(*args, **kwargs)

    def after(self, *args, **kwargs):
        self.window.after(*args, **kwargs)


# single widgets


class HOPSLabel(HOPSWidget):

    def __init__(self, window, text=' ', instance=str):

        if instance == str:
            self.variable = StringVar(window.root, value=text)
        elif instance == float:
            self.variable = DoubleVar(window.root, value=text)
        elif instance == int:
            self.variable = IntVar(window.root, value=text)
        elif instance == bool:
            self.variable = BooleanVar(window.root, value=text)

        widget = Label(window.root, textvar=self.variable)

        HOPSWidget.__init__(self, window, widget, 'Label')

    def set(self, value):
        self.variable.set(value)

    def get(self):
        return self.variable.get()


class HOPSEntry(HOPSWidget):

    def __init__(self, window, value=' ', instance=str, command=None):

        self.instance = instance

        self.variable = StringVar(window.root, value=value)

        widget = Entry(window.main_frame, textvar=self.variable, validate='key')

        widget['validatecommand'] = (widget.register(self.validate), '%P', '%d')

        if command:
            widget.bind(sequence='<KeyRelease>', func=command)

        HOPSWidget.__init__(self, window, widget, 'Entry')

    def get(self):

        if self.instance == str:
            return self.variable.get()
        else:
            if self.variable.get() == '':
                return 0
            else:
                return self.instance(self.variable.get())

    def validate(self, input_str, typing):

        if typing == '1':
            try:
                _ = self.instance(input_str)
                return True
            except ValueError:
                return False

        else:
            return True

    def set(self, value):
        self.variable.set(value)



class HOPSButton(HOPSWidget):

    def __init__(self, window, text=' ', command=None, **kwargs):

        self.variable = StringVar(window.root, value=text)
        self.command=command

        widget = Button(window.main_frame, textvar=self.variable, command=self.internal_command, **kwargs)

        HOPSWidget.__init__(self, window, widget, 'Button')

    def internal_command(self):

        if isinstance(self.command, list):
            functions = [self.nested_after(self.command[-2], self.command[-1])]
            for i in range(len(self.command)-2):

                functions.append(self.nested_after(self.command[-i-3], functions[-1]))
            self.window.after(functions[-1])

        else:
            self.command()

    def nested_after(self, function1, function2):

        def xx():
            function1()
            self.window.after(function2, 10)

        return xx


class HOPSCheckButton(HOPSWidget):

    def __init__(self, window, text=' ', initial=False, command=None):

        self.variable = BooleanVar(window.root, value=initial)

        if command:
            widget = Checkbutton(window.main_frame, text=text, variable=self.variable, command=command)
        else:
            widget = Checkbutton(window.main_frame, text=text, variable=self.variable)

        HOPSWidget.__init__(self, window, widget, 'Checkbutton')

    def set(self, value):
        self.variable.set(value)

    def get(self):
        return self.variable.get()


class HOPSDropDown(HOPSWidget):

    def __init__(self, window, initial, options, instance, command=None, width=None):

        try:
            initial = instance(initial)
            if initial not in options:
                initial = options[0]
        except:
            initial = options[0]

        self.variable = StringVar(window.root, value=initial)
        self.instance = instance
        self.command = command

        try:

            combostyle = Style()
            combostyle.theme_create('combostyle', parent='alt',
                                    settings={'TCombobox': {'configure':
                                                                {'selectbackground': 'white',
                                                                 'fieldbackground': 'white',
                                                                 'background': 'white'}}})
            combostyle.theme_use('combostyle')

        except:
            pass
        if width:
            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=width)
        else:
            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=int(window.log.entries_width * 0.8))

        options = [str(ff) for ff in options]
        widget['values'] = tuple(options)

        widget.bind('<<ComboboxSelected>>', self.internal_command)

        HOPSWidget.__init__(self, window, widget, 'DropDown')

    def internal_command(self, *entry):
        self.widget.selection_clear()
        if self.command:
            self.command()

    def set(self, value):
        self.variable.set(self.instance(value))

    def get(self):
        return self.instance(self.variable.get())


# multi-widgets


class HOPSListDisplay(HOPSWidget):

    def __init__(self, window, font='Courier'):

        widget = Frame(window.main_frame)

        scrollbar = Scrollbar(widget)
        scrollbar.pack(side=RIGHT, fill=Y)
        self.listbox = Listbox(widget, yscrollcommand=scrollbar.set, font=font)
        self.listbox.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.config(command=self.listbox.yview)

        window.root.columnconfigure(0, weight=1)
        window.root.rowconfigure(0, weight=1)

        HOPSWidget.__init__(self, window, widget, 'ListDisplay')

    def update_list(self, list_to_add):

        self.listbox.delete(0, END)
        for row in list_to_add:
            self.listbox.insert(END, str(row))

    def grid(self, *args, **kwargs):
        return self.widget.grid(*args, **kwargs, sticky=E+W+N+S)


class HOPSProgressbar(HOPSWidget):

    def __init__(self, window, task='Process', length=0.2):

        self.task = task

        widget = Frame(window.main_frame)
        self.progress = DoubleVar(widget, value=0)
        self.progressbar = Progressbar(widget, variable=self.progress,
                                       length=length * window.root.winfo_screenwidth(), orient=HORIZONTAL,
                                       maximum=100, mode='determinate', value=0)

        self.task_label = Label(widget, text=self.task + ' - 000.0 % - time left: 00:00:00')
        self.element_label = Label(widget, text=' ')

        self.task_label.pack()
        self.progressbar.pack()
        self.element_label.pack()

        self.start_time = 0
        self.skip_time = 0
        self.total_iterations = 1
        self.current_iteration = 0
        self.show_at = [0]

        HOPSWidget.__init__(self, window, widget, 'Progressbar')

    def reset(self):

        self.progress.set(0)
        self.task_label.configure(text=self.task + ' - 000.0 % - time left: 00:00:00')

    def initiate(self, total_iterations, show_every=1):

        if isinstance(total_iterations, list):
            self.total_iterations = len(total_iterations)
            self.elements = total_iterations
        else:
            self.total_iterations = int(total_iterations)
            self.elements = None

        self.current_iteration = 0
        self.start_time = time.time()

        self.show_at = list(range(1, self.total_iterations, int(show_every))) + [int(self.total_iterations)]

    def show_message(self, message):
        if message != self.element_label['text']:
            self.element_label.configure(text=message)

    def update(self, step=1, skip=0):

        self.current_iteration += step
        self.skip_time += skip

        delta_time = time.time() - self.start_time - self.skip_time

        time_left = str(datetime.timedelta(
            seconds=int((self.total_iterations - self.current_iteration) * delta_time / self.current_iteration)))

        percent = round(100 * float(self.current_iteration) / float(self.total_iterations), 1)

        if self.current_iteration in self.show_at :
            self.task_label.configure(text=self.task + ' - {0} % - time left: {1} '.format(percent, time_left))
            if self.elements:
                self.element_label.configure(text=self.elements[self.current_iteration - 1])

            self.progress.set(percent)


class HOPSFigureWindow(HOPSWidget):

    def __init__(self, window,
                 figsize=None, max_figsize_percent=(0.9,0.8),
                 show_nav=False, subplots_adjust=None,):

        widget = Frame(window.main_frame)

        self.max_figsize_percent = max_figsize_percent

        if figsize:
            self.figure = matplotlib.figure.Figure(figsize=figsize)
        else:
            self.figure = matplotlib.figure.Figure()

        self.figure.patch.set_facecolor('white')
        self.canvas = FigureCanvasTkAgg(self.figure, master=widget)
        self.canvas.get_tk_widget().pack(side=TOP)

        if subplots_adjust:
            self.figure.subplots_adjust(left=subplots_adjust[0], right=subplots_adjust[1],
                                        bottom=subplots_adjust[2], top=subplots_adjust[3])

        if show_nav:
            toolbar = NavigationToolbar2Tk(self.canvas, widget)
            toolbar.pack(side=BOTTOM)

        HOPSWidget.__init__(self, window, widget, 'FigureWindow')

    def adjust_size(self):

        self.window.reposition()

        current_figure_height = self.canvas.get_tk_widget().winfo_reqheight()
        current_figure_width = self.canvas.get_tk_widget().winfo_width()

        if self.root.winfo_reqheight() > self.max_figsize_percent[1] * self.root.winfo_screenheight():
            dh = (self.root.winfo_reqheight() - self.max_figsize_percent[1] * self.root.winfo_screenheight())
            new_figure_height = max(current_figure_height - dh, 100)
            self.canvas.get_tk_widget().config(height=new_figure_height)

        if self.root.winfo_reqwidth() > self.max_figsize_percent[0] * self.root.winfo_screenwidth():
            dw = (self.root.winfo_reqwidth() - self.max_figsize_percent[0] * self.root.winfo_screenwidth())
            new_figure_width = max(current_figure_width - dw, 100)
            self.canvas.get_tk_widget().config(width=new_figure_width)

        self.window.reposition()

    def draw(self, update_level=1):
        self.canvas.draw()
        if update_level == 0:
            pass
        elif update_level == 1:
            self.window.update_idletasks()
        elif update_level == 2:
            self.window.update()


class HOPSFitsWindow(HOPSWidget):

    def __init__(self, window, fits_data=None, fits_header=None,
                 input_name=None, input_options=None,
                 max_figsize_percent=(0.9,0.8),  subplots_adjust=None,
                 show_nav=False, show_controls=False, show_axes=False):

        widget = Frame(window.root)
        HOPSWidget.__init__(self, window, widget, 'FitsWindow')

        self.show_axes = show_axes
        self.show_half = False
        self.max_figsize_percent = max_figsize_percent

        self.figure = matplotlib.figure.Figure()
        self.figure.patch.set_facecolor('white')
        self.canvas = FigureCanvasTkAgg(self.figure, widget)
        self.ax = self.figure.add_subplot(111)
        self.ax.tick_params(axis='y', rotation=90)

        if not self.show_axes:
            self.ax.axis('off')

        if subplots_adjust:
            self.figure.subplots_adjust(left=subplots_adjust[0], right=subplots_adjust[1],
                                        bottom=subplots_adjust[2], top=subplots_adjust[3])
        else:
            self.figure.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

        self.data = None
        self.mean = 0
        self.std = 0
        self.image = None
        self.sqrt_vmin = DoubleVar(widget, value=0)
        self.vmin = IntVar(widget, value=0)
        self.sqrt_vmax = DoubleVar(widget, value=10000)
        self.vmax = IntVar(widget, value=10000)
        self.gamma = DoubleVar(widget, value=0)
        self.flip = IntVar(widget, value=0)
        self.mirror = IntVar(widget, value=0)
        self.white_sky = IntVar(widget, value=0)

        self.fits_name = StringVar(widget, value=input_name)
        self.fits_name_label = Label(widget, textvar=self.fits_name)

        # extra widgets

        control_frame = Frame(widget)
        self.control_frame = control_frame

        self.info_label = Label(control_frame, text='Scroll up/down to zoom in/out. Click & drag to move the image.')

        self.mouse_data = StringVar(control_frame, value=' ')
        self.mouse_data_label = Label(control_frame, textvar=self.mouse_data)

        self.black_entry = Scale(control_frame, resolution=0.1, variable=self.sqrt_vmin, orient=HORIZONTAL, showvalue=False)
        self.black_entry.bind("<B1-Motion>", self.contrast)
        self.black_entry.bind("<ButtonRelease-1>", self.contrast)
        self.black_entry_label_0 = Label(control_frame, text='Minimum = ', anchor=E)
        self.black_entry_label = Label(control_frame, textvar=self.vmin, anchor=W)
        self.black_entry['from_'] = 1
        self.black_entry['to'] = 1000

        self.white_entry = Scale(control_frame, resolution=0.1, variable=self.sqrt_vmax, orient=HORIZONTAL, showvalue=False)
        self.white_entry.bind("<B1-Motion>", self.contrast)
        self.white_entry.bind("<ButtonRelease-1>", self.contrast)
        self.white_entry_label_0 = Label(control_frame, text='Maximum = ', anchor=E)
        self.white_entry_label = Label(control_frame, textvar=self.vmax, anchor=W)
        self.white_entry['from_'] = 1
        self.white_entry['to'] = 1000

        self.gamma_entry = Scale(control_frame, resolution=0.001, variable=self.gamma, orient=HORIZONTAL, showvalue=False)
        self.gamma_entry.bind("<B1-Motion>", self.contrast)
        self.gamma_entry.bind("<ButtonRelease-1>", self.contrast)
        self.gamma_entry_label_0 = Label(control_frame, text='Stretch factor = ', anchor=E)
        self.gamma_entry_label = Label(control_frame, textvar=self.gamma, anchor=W)
        self.gamma_entry['from_'] = 0
        self.gamma_entry['to'] = 1

        self.flip_button = Checkbutton(control_frame, text='Flip', variable=self.flip, command=self.flip_fov)
        self.mirror_button = Checkbutton(control_frame, text='Mirror', variable=self.mirror, command=self.mirror_fov)
        self.reverse_color_button = Checkbutton(control_frame, text='White Sky', variable=self.white_sky, command=self.reverse_color)
        self.reset_button = Button(control_frame, text='RESET', command=self.reset)

        self.info_label.grid(row=1, column=1, columnspan=4)
        self.mouse_data_label.grid(row=2, column=1, columnspan=4)
        self.black_entry_label_0.grid(row=3, column=1, columnspan=2)
        self.black_entry_label.grid(row=3, column=3)
        self.black_entry.grid(row=4, column=1, columnspan=4, sticky=N+S+E+W)
        self.white_entry_label_0.grid(row=5, column=1, columnspan=2)
        self.white_entry_label.grid(row=5, column=3)
        self.white_entry.grid(row=6, column=1, columnspan=4, sticky=N+S+E+W)
        self.gamma_entry_label_0.grid(row=7, column=1, columnspan=2)
        self.gamma_entry_label.grid(row=7, column=3)
        self.gamma_entry.grid(row=8, column=1, columnspan=4, sticky=N+S+E+W)
        self.reset_button.grid(row=9, column=1)
        self.flip_button.grid(row=9, column=2)
        self.mirror_button.grid(row=9, column=3)
        self.reverse_color_button.grid(row=9, column=4)
        Label(control_frame, text=' ').grid(row=10, column=1, columnspan=4)

        self.picked = False

        self.canvas.get_tk_widget().pack(side=TOP, fill='both', expand=True)
        if show_nav:
            toolbar = NavigationToolbar2Tk(self.canvas, self.widget)
            toolbar.pack(side=BOTTOM)

        self.fits_name_label.pack()

        if show_controls:
            control_frame.pack()
            self.canvas.callbacks.connect('scroll_event', self.zoom)
            self.canvas.callbacks.connect('motion_notify_event', self.move)
            self.canvas.callbacks.connect('button_press_event', self.pick)
            self.canvas.callbacks.connect('button_release_event', self.pick)
        else:
            self.canvas.stop_event_loop

        if type(fits_data) != type(None):
            self.load_fits(fits_data, fits_header, input_name, input_options)

    def adjust_size(self):

        self.window.reposition()

        current_figure_height = self.canvas.get_tk_widget().winfo_reqheight()
        current_figure_width = self.canvas.get_tk_widget().winfo_width()

        if self.root.winfo_reqheight() > self.max_figsize_percent[1] * self.root.winfo_screenheight():
            dh = (self.root.winfo_reqheight() - self.max_figsize_percent[1] * self.root.winfo_screenheight())
            new_figure_height = max(current_figure_height - dh, 100)
            self.canvas.get_tk_widget().config(height=new_figure_height)

        if self.root.winfo_reqwidth() > self.max_figsize_percent[0] * self.root.winfo_screenwidth():
            dw = (self.root.winfo_reqwidth() - self.max_figsize_percent[0] * self.root.winfo_screenwidth())
            new_figure_width = max(current_figure_width - dw, 100)
            self.canvas.get_tk_widget().config(width=new_figure_width)

        self.window.reposition()

    def load_fits(self, fits_data, fits_header, input_name=None, input_options=None, draw=True, shift=0, show_half=False):

        self.show_half = show_half

        if input_name:
            input_name = os.path.split(input_name)[1]
            if len(input_name) > 50:
                split = [input_name[i:i + 50] for i in range(0, len(input_name), 50)]
                input_name = '\n'.join(split)

        self.fits_name.set(input_name)

        self.data = np.ones_like(fits_data) * fits_data

        try:
            self.mean = fits_header[self.window.log.mean_key]
            self.std = fits_header[self.window.log.std_key]
        except:
            self.mean, self.std = plc.mean_std_from_median_mad(fits_data)

        self.black_entry['from_'] = np.sqrt(max(0, np.min(self.data)))
        self.black_entry['to'] = np.sqrt(np.max(self.data))

        self.white_entry['from_'] = np.sqrt(max(0, np.min(self.data)))
        self.white_entry['to'] = np.sqrt(np.max(self.data))

        self.ax.cla()
        if not self.show_axes:
            self.ax.axis('off')
        self.ax.tick_params(axis='y', rotation=90)

        self.vmin.set(max(1, int(self.mean + self.window.log.frame_low_std * self.std)))
        self.vmax.set(max(1, int(self.mean + self.window.log.frame_upper_std * self.std)))
        self.gamma.set(0)

        if input_options:
            if input_options[0] != 'auto':
                self.vmin.set(max(1, int(self.mean + input_options[0] * self.std)))
            if input_options[1] != 'auto':
                self.vmax.set(max(1, int(self.mean + input_options[1] * self.std)))
            self.gamma.set(input_options[2])
            self.flip.set(input_options[3])
            self.mirror.set(input_options[4])
            self.white_sky.set(input_options[5])

        self.sqrt_vmin.set(np.sqrt(self.vmin.get()))
        self.sqrt_vmax.set(np.sqrt(self.vmax.get()))

        if self.show_half:
            xl = len(self.data[0])
            yl = len(self.data)
            x1, x2 = int(0.25*xl), int(0.75*xl)
            y1, y2 = int(0.25*yl), int(0.75*yl)
            self.image = self.ax.imshow(self.data[y1:y2,x1:x2] ** (10 ** - self.gamma.get()), origin='lower',
                                        extent=(x1, x2, y1, y2),
                                        cmap=Greys, vmin=self.vmin.get(), vmax=self.vmax.get())
        else:
            self.image = self.ax.imshow(self.data ** (10 ** - self.gamma.get()), origin='lower',
                                        extent=(0, len(self.data[0]), 0, len(self.data)),
                                        cmap=Greys, vmin=self.vmin.get(), vmax=self.vmax.get())

        if self.white_sky.get():
            self.image.set_cmap(Greys)
        else:
            self.image.set_cmap(Greys_r)

        if self.flip.get():
            self.ax.set_ylim(len(self.data), 0)
        else:
            self.ax.set_ylim(0, len(self.data))

        if self.mirror.get():
            self.ax.set_xlim(len(self.data[0]), 0)
        else:
            self.ax.set_xlim(0, len(self.data[0]))

        try:
            if 'auto' not in input_options[6:]:
                self.ax.set_xlim(*input_options[6:8])
                self.ax.set_ylim(*input_options[8:10])
        except Exception as e:
            pass

        if draw:
            self.draw()

    def reverse_color(self):
        if self.white_sky.get():
            self.image.set_cmap(Greys)
        else:
            self.image.set_cmap(Greys_r)

        self.canvas.draw()

    def flip_fov(self):
        lims = self.ax.get_ylim()
        if self.flip.get():
            self.ax.set_ylim(max(lims), min(lims))
        else:
            self.ax.set_ylim(min(lims), max(lims))
        self.canvas.draw()

    def mirror_fov(self):
        lims = self.ax.get_xlim()
        if self.mirror.get():
            self.ax.set_xlim(max(lims), min(lims))
        else:
            self.ax.set_xlim(min(lims), max(lims))
        self.canvas.draw()

    def contrast(self, event):

        if self.sqrt_vmin.get() >= self.sqrt_vmax.get():
            self.sqrt_vmin.set(self.sqrt_vmax.get() - 1)

        self.vmin.set(int(self.sqrt_vmin.get() ** 2))
        self.vmax.set(int(self.sqrt_vmax.get() ** 2))

        self.image.set_data(np.maximum(0, self.data) ** (10 ** - self.gamma.get()))

        self.image.set_clim(self.vmin.get() ** (10 ** -self.gamma.get()), self.vmax.get() ** (10 ** -self.gamma.get()))
        self.canvas.draw()

    def get_fov_options(self):
        try:
            return [(self.vmin.get() - self.mean)/self.std, (self.vmax.get() - self.mean)/self.std,
                    self.gamma.get(), self.flip.get(), self.mirror.get(), self.white_sky.get(),
                    int(self.ax.get_xlim()[0]), int(self.ax.get_xlim()[1]),
                    int(self.ax.get_ylim()[0]), int(self.ax.get_ylim()[1])]
        except:
            return None

    def pick(self, event):

        if isinstance(event, matplotlib.backend_bases.MouseEvent):

            if event.inaxes is None:
                pass

            elif event.dblclick:
                pass

            elif event.name == 'button_press_event':

                self.picked = (event.xdata, event.ydata)

            elif event.name == 'button_release_event':

                self.picked = False

    def move(self, event):

        if isinstance(event, matplotlib.backend_bases.MouseEvent):

            if event.inaxes is None:
                pass

            elif event.name == 'motion_notify_event':

                try:
                    self.mouse_data.set('Mouse on: x={0:.2f}, y={1:.2f}, counts={2:.2f}'.format(
                        event.xdata, event.ydata, self.data[int(event.ydata), int(event.xdata)]))
                except:
                    self.mouse_data.set('Mouse on: x={0:.2f}, y={1:.2f}, counts={2}'.format(
                        event.xdata, event.ydata, '-'))

                if self.picked:

                    dx = event.xdata - self.picked[0]
                    dy = event.ydata - self.picked[1]

                    self.ax.set_xlim(self.ax.get_xlim()[0] - dx, self.ax.get_xlim()[1] - dx)
                    self.ax.set_ylim(self.ax.get_ylim()[0] - dy, self.ax.get_ylim()[1] - dy)

                    self.canvas.draw()

    def zoom(self, event):

        if isinstance(event, matplotlib.backend_bases.MouseEvent):

            if event.inaxes is None:
                pass

            elif event.name == 'scroll_event':

                zoom_factor = 1.2
                scale_factor = 1.0

                if event.button == 'up':
                    scale_factor = 1 / zoom_factor
                elif event.button == 'down':
                    scale_factor = zoom_factor

                xdata = event.xdata
                ydata = event.ydata

                cur_xlim = self.ax.get_xlim()
                cur_ylim = self.ax.get_ylim()

                cur_xrange = (cur_xlim[1] - cur_xlim[0])
                cur_yrange = (cur_ylim[1] - cur_ylim[0])

                new_xrange = cur_xrange * scale_factor
                new_yrange = cur_yrange * scale_factor

                new_xmin = xdata - new_xrange * (xdata - cur_xlim[0]) / cur_xrange
                new_ymin = ydata - new_yrange * (ydata - cur_ylim[0]) / cur_yrange

                self.ax.set_xlim([new_xmin, new_xmin + new_xrange])
                self.ax.set_ylim([new_ymin, new_ymin + new_yrange])

                self.canvas.draw()

    def reset(self):

        self.ax.set_xlim(0, len(self.data[0]) + 5)
        self.ax.set_ylim(0, len(self.data) + 5)

        if self.flip.get():
            self.ax.set_ylim(self.ax.get_ylim()[1], self.ax.get_ylim()[0])

        if self.mirror.get():
            self.ax.set_xlim(self.ax.get_xlim()[1], self.ax.get_xlim()[0])

        self.vmin.set(max(1, int(self.mean + self.window.log.frame_low_std * self.std)))
        self.sqrt_vmin.set(np.sqrt(self.vmin.get()))
        self.vmax.set(max(1, int(self.mean + self.window.log.frame_upper_std * self.std)))
        self.sqrt_vmax.set(np.sqrt(self.vmax.get()))

        self.gamma.set(0)

        self.image.set_data(np.maximum(0, self.data) ** (10 ** -self.gamma.get()))

        self.image.set_clim(self.vmin.get() ** (10 ** -self.gamma.get()), self.vmax.get() ** (10 ** -self.gamma.get()))
        self.draw()

    def draw(self, update_level=1):
        self.canvas.draw()

    def disable(self):
        for child in self.widget.winfo_children():
            try:
                child.configure(state='disable')
            except:
                pass

    def activate(self):
        for child in self.widget.winfo_children():
            try:
                child.configure(state='active')
            except:
                pass
