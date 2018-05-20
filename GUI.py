import tkinter as tk
from tkinter import ttk
import matplotlib
import os
from plotting_class import FilePlotting

LARGE_FONT = ("Verdana", 10)


class MDAnalysis(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Analysis of MD simulations")
        container = ttk.Frame(self, width=200, height=200)

        self.frames = {}

        for F in (MainPage, PageOne):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(MainPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.particles = tk.IntVar()
        self.steps = tk.IntVar()
        # plotting entries
        energies_label = ttk.Label(text="Energies", font=LARGE_FONT)
        energies_label.grid(row=0, column=0, columnspan=2)

        stat_label = ttk.Label(text="Statistical Analysis", font=LARGE_FONT)
        stat_label.grid(row=0, column=3, columnspan=3)

        # Entry fields
        # self. included in these 2 variables because otherwise it was not working
        self.steps.set(10000)
        __steps_label = ttk.Label(text="N(steps)")
        __steps_label.grid(row=0, column=6)
        __steps = ttk.Entry(textvariable=self.steps, width=6)
        __steps.grid(row=0, column=7)

        self.particles.set(1000)
        __particles_label = ttk.Label(text="Particles")
        __particles_label.grid(row=0, column=8)
        __particles = ttk.Entry(textvariable=self.particles, width=4)
        __particles.grid(row=0, column=9)

        rho = tk.DoubleVar()
        __rho_label = ttk.Label(text="\u03C1")
        __rho_label.grid(row=0, column=10)
        __rho = ttk.Entry(textvariable=rho, width=6)
        __rho.grid(row=0, column=11)

        t = tk.DoubleVar()
        __t_label = ttk.Label(text="T")
        __t_label.grid(row=0, column=12)
        __t = ttk.Entry(textvariable=t, width=6)
        __t.grid(row=0, column=13)

        a = tk.DoubleVar()
        __a_label = ttk.Label(text="a")
        __a_label.grid(row=0, column=14)
        __a = ttk.Entry(textvariable=a, width=6)
        __a.grid(row=0, column=15)

        n = tk.IntVar()
        __n_label = ttk.Label(text="n")
        __n_label.grid(row=0, column=16)
        __n = ttk.Entry(textvariable=n, width=2)
        __n.grid(row=0, column=17)

        obj = FilePlotting(self.steps.get(), self.particles.get())
        os.chdir('../../Archives of Data/')

        # Buttons
        # Energies vs time (for single run
        u_en = ttk.Button(text="U")
        u_en.grid(row=1, column=0)

        k_en = ttk.Button(text="K")
        k_en.grid(row=1, column=1)

        total_en = ttk.Button(text="U+K")
        total_en.grid(row=2, column=0)

        all_en = ttk.Button(text="All")
        all_en.grid(row=2, column=1)

        # Statistical Quantities
        # tkinter variables not converted to Python variables

        button_rdf = ttk.Button(text="RDF",
                                command=lambda: obj.rdf(rho.get(), t.get(), n.get(), a.get()))
        button_rdf.grid(row=1, column=3)

        msd = ttk.Button(text="MSD",
                         command=lambda: obj.msd(rho.get(), t.get(), n.get(), a.get()))
        msd.grid(row=1, column=4)

        vaf = ttk.Button(text="VAF",
                         command=lambda: obj.vaf(rho.get(), t.get(), n.get(), a.get()))
        vaf.grid(row=1, column=5)


class PageOne(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Page 1", font=LARGE_FONT)
        label.grid(row=0, column=0, pady=10, padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(MainPage))
        button1.grid(row=1, column=0)


app = MDAnalysis()
app.mainloop()
