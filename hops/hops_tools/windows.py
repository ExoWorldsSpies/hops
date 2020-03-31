
import os
from tkinter import Tk, TclError
from tkinter import Label, Button, Entry, Checkbutton, Scrollbar, Listbox, PhotoImage, Radiobutton, Scale, Frame
from tkinter import StringVar, BooleanVar, DoubleVar, IntVar
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, BOTH, Y, HORIZONTAL

from tkinter.ttk import Combobox, Style, Progressbar
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showinfo

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    NavigationToolbar2Tk = NavigationToolbar2TkAgg
# NavigationToolbar2TkAgg

import webbrowser

from hops.hops_tools.logs import log


def openweb():
    webbrowser.open("https://www.exoworldsspies.com/en/software", new=1)

def openweb_simbad(radec_string):

    def mock(radec_string=radec_string):
        radec_string = radec_string.replace('+', '%2B').replace(' ', '+')
        webbrowser.open("http://simbad.u-strasbg.fr/simbad/sim-coo?Coord={0}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=20&Radius.unit=arcmin&submit=submit+query&CoordList=".format(radec_string), new=1)

    return mock


class HOPSWindow:

    def __init__(self, name, sizex=None, sizey=None, position=5):

        self.root = Tk()
        self.root.wm_title(name)
        if sizex and sizey:
            self.root.geometry('{0}x{1}'.format(int(self.root.winfo_screenwidth() / sizex),
                                                int(self.root.winfo_screenheight() / sizey)))

        self.hide()

        self.finalised = False
        self.position = position

        self.running = self.BooleanVar(False)
        self.exit = False

    def exit_command(self):
        self.exit = True

    def exit_python(self):
        os._exit(-1)

    def no_action(self):
        pass

    def finalise(self):

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

        self.root.lift()
        self.root.wm_attributes("-topmost", 1)
        self.root.after_idle(self.root.attributes, '-topmost', 0)

    def hide(self):

        self.root.withdraw()

    def show(self):

        if not self.finalised:
            self.finalise()
            self.finalised = True

        self.root.deiconify()

    def loop(self):

        if not self.finalised:
            self.finalise()
            self.finalised = True

        self.root.deiconify()

        self.root.mainloop()

    def close(self):

        self.root.quit()
        self.root.destroy()

    def setup_window(self, objects, title_font='times', main_font='times', button_font='times', entries_wd=20, entries_bd=3, buttons_bd=5):

        # print(font.families())
        font_size = min(15, self.root.winfo_screenheight()/60)

        button_font = (button_font, int(1.1 * font_size), 'bold')
        main_font = (main_font, int(font_size))
        title_font = (title_font, int(1.5 * font_size), 'bold')

        for row in range(len(objects)):
            if len(objects[row]) == 0:
                label_empty = Label(self.root, text='')
                label_empty.grid(row=row, column=100)
            else:
                for obj in objects[row]:

                    if obj[0].winfo_class() == 'Button':
                        obj[0].config(borderwidth=buttons_bd, font=button_font, padx=1, pady=1)
                    elif obj[0].winfo_class() == 'Entry':
                        obj[0].configure(width=entries_wd, bd=entries_bd, font=main_font)
                    elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
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

    def update(self):

        self.root.update()

    def update_idletasks(self):

        self.root.update_idletasks()

    def after(self, *args):

        self.root.after(*args)

    def Label(self, *args, **kwargs):

        return Label(self.root, *args, **kwargs)

    def Entry(self, *args, **kwargs):

        return Entry(self.root, *args, **kwargs)

    def Button(self, *args, **kwargs):

        return Button(self.root, *args, **kwargs)

    def Checkbutton(self, *args, **kwargs):

        return Checkbutton(self.root, *args, **kwargs)

    def Scrollbar(self, *args, **kwargs):

        return Scrollbar(self.root, *args, **kwargs)

    def Listbox(self, *args, **kwargs):

        return Listbox(self.root, *args, **kwargs)

    def Radiobutton(self, *args, **kwargs):

        return Radiobutton(self.root, *args, **kwargs)

    def Scale(self, *args, **kwargs):

        return Scale(self.root, *args, **kwargs)

    def Combobox(self, *args, **kwargs):

        return Combobox(self.root, *args, **kwargs)

    def Style(self):

        return Style()

    def Progressbar(self):

        return Progressbar(self.root, orient=HORIZONTAL, length=self.root.winfo_screenwidth() / 5.0,
                           maximum=100, mode='determinate', value=0)

    def Frame(self, *args, **kwargs):

        return Frame(self.root, *args, **kwargs)

    def FigureCanvasTkAgg(self, figure):
        return FigureCanvasTkAgg(figure, self.root)

    def NavigationToolbar2Tk(self, canvas):
        return NavigationToolbar2Tk(canvas, self.root)

    def StringVar(self, value):

        return StringVar(self.root, value=value)

    def BooleanVar(self, value):

        return BooleanVar(self.root, value=value)

    def DoubleVar(self, value):

        return DoubleVar(self.root, value=value)

    def IntVar(self, value):

        return IntVar(self.root, value=value)

    def askdirectory(self):

        return askdirectory()


class MainWindow(HOPSWindow):

    def __init__(self, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self, name, sizex, sizey, position)
        self.root.protocol('WM_DELETE_WINDOW', self.exit_python)

        self.close_buttons = []
        for i in range(20):
            photo2 = PhotoImage(file=log.close_icon)
            self.close_buttons.append(Button(image=photo2))
            # this is important for photoimage, needs a reference
            self.close_buttons[-1].image = photo2

        self.created_by_label = self.Label(text=log.read_log('windows', 'created_by').replace(',', '\n'))
        self.Btn = self.Button(text="HOPS UPDATES &\nUSER MANUAL", command=openweb)

        self.observing_planner_button = self.Button(text='OBSERVATION\nPLANNER')

        self.my_profile_window = ProfileWindow()
        self.my_profile_button = self.Button(text='MY PROFILE', command=self.my_profile_window.show)

    def get_logo(self):
        photo = PhotoImage(file=log.holomon_logo)
        logo = Label(image=photo)
        # this is important for photoimage, needs a reference
        logo.image = photo
        # del photo
        return logo


class SideWindow(HOPSWindow):

    def __init__(self, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self, name, sizex, sizey, position)
        self.root.protocol('WM_DELETE_WINDOW', self.no_action)


class AddOnWindow(HOPSWindow):

    def __init__(self, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self,name, sizex, sizey, position)
        self.root.protocol('WM_DELETE_WINDOW', self.hide)


class ProgressWindow(HOPSWindow):

    def __init__(self, name, sizex=None, sizey=None, position=5):

        HOPSWindow.__init__(self,name, sizex, sizey, position)
        self.root.protocol('WM_DELETE_WINDOW', self.exit_command)


class ProfileWindow(AddOnWindow):

    def __init__(self):

        AddOnWindow.__init__(self, 'My Profile', 2, 1.1, 1)

        core_headers = {ff.split(':')[0]: log.read_log_profile(ff.split(':')[0])
                        for ff in open(log.files['log_profile'], 'r').readlines()}

        local_headers = log.open_yaml(log.files['local_log_profile'])

        variables = {}
        labels = {}
        entries = {}
        for row, header in enumerate(core_headers):
            if header in local_headers:
                variables[header] = self.StringVar(local_headers[header])
                labels[header] = self.Label(text=header)
                entries[header] = self.Entry(textvariable=variables[header])
            else:
                variables[header] = self.StringVar(core_headers[header])
                labels[header] = self.Label(text=header)
                entries[header] = self.Entry(textvariable=variables[header])

        def update_headers():
            new_local_profile = {}
            for header2 in variables:
                new_local_profile[header2] = variables[header2].get()
            log.save_yaml(new_local_profile, log.files['local_log_profile'])

        update_headers_button = self.Button(text='UPDATE')
        update_headers_button['command'] = update_headers

        stucture = [[], [[update_headers_button, 2]], []]
        for header in list(core_headers.keys())[:int(len(list(core_headers.keys())) / 2)]:
            stucture.append([[labels[header], 1], [entries[header], 2]])

        stucture.append([])

        for jj, header in enumerate(list(core_headers.keys())[int(len(list(core_headers.keys())) / 2):]):
            stucture[3 + jj].append([labels[header], 3])
            stucture[3 + jj].append([entries[header], 4])

        stucture.append([])

        self.setup_window(stucture)


def setup_window(window, objects, title_font=None, main_font=None, button_font=None, entries_bd=3, buttons_bd=5):
    screenheigth = window.winfo_screenheight()

    if button_font is None:
        button_font = ['times', int(screenheigth/55), 'bold']

    if main_font is None:
        main_font = ['times', int(screenheigth/60)]

    if title_font is None:
        title_font = ['times', int(screenheigth/40), 'bold']

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].config(borderwidth=buttons_bd, font=button_font, padx=3, pady=3)
                elif obj[0].winfo_class() == 'Entry':
                    obj[0].configure(bd=entries_bd, font=main_font)
                elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
                    if len(obj) == 5:
                        if obj[4] == 'title':
                            obj[0].configure(font=title_font)
                        else:
                            obj[0].configure(font=main_font)
                    else:
                        obj[0].configure(font=main_font)

                if len(obj) >= 4:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                elif len(obj) == 3:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                else:
                    obj[0].grid(row=row, column=obj[1])
