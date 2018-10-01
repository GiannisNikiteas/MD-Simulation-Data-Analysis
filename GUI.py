import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
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

        # Entry fields
        # self. included in these 2 variables because otherwise it was not working
        self.steps.set(10000)
        __steps_label = ttk.Label(text="N(steps)")
        __steps_label.grid(row=0, column=6)
        __steps = ttk.Entry(textvariable=self.steps, width=6)
        __steps.grid(row=0, column=7, pady=5)

        self.particles.set(1000)
        __particles_label = ttk.Label(text="Particles")
        __particles_label.grid(row=0, column=8)
        __particles = ttk.Entry(textvariable=self.particles, width=4)
        __particles.grid(row=0, column=9, pady=5)

        rho = tk.DoubleVar()
        rho.set(0.5)
        __rho_label = ttk.Label(text="\u03C1")
        __rho_label.grid(row=0, column=10)
        __rho = ttk.Entry(textvariable=rho, width=6)
        __rho.grid(row=0, column=11, pady=5)

        t = tk.DoubleVar()
        t.set(0.5)
        __t_label = ttk.Label(text="T")
        __t_label.grid(row=0, column=12)
        __t = ttk.Entry(textvariable=t, width=6)
        __t.grid(row=0, column=13, pady=5)

        a = tk.DoubleVar()
        a.set(0.5)
        __a_label = ttk.Label(text="a")
        __a_label.grid(row=0, column=14)
        __a = ttk.Entry(textvariable=a, width=6)
        __a.grid(row=0, column=15, pady=5)

        n = tk.IntVar()
        n.set(12)
        __n_label = ttk.Label(text="n")
        __n_label.grid(row=0, column=16)
        __n = ttk.Entry(textvariable=n, width=2)
        __n.grid(row=0, column=17, pady=5)

        obj = FilePlotting(self.steps.get(), self.particles.get())
        os.chdir('../../Archives of Data/')

        # Buttons
        # plotting entries
        energies_label = ttk.Label(text="Energies", font=LARGE_FONT)
        energies_label.grid(row=0, column=0, columnspan=2, pady=5)

        stat_label = ttk.Label(text="Statistical Analysis", font=LARGE_FONT)
        stat_label.grid(row=0, column=3, columnspan=3, pady=5)

        # todo: this is a hack. find better sol'n
        empty_label = ttk.Label(text=" ")
        empty_label.grid(row=0, column=2)
        self.grid_columnconfigure(2, minsize=1000)

        # Energies vs time for single runs
        u_en = ttk.Button(text="U",
                          command=obj.potential_data(rho.get(), t.get(), n.get(), a.get()))
        u_en.grid(row=1, column=0)

        k_en = ttk.Button(text="K")
        k_en.grid(row=1, column=1)

        total_en = ttk.Button(text="U+K")
        total_en.grid(row=2, column=0)

        all_en = ttk.Button(text="All",
                            command=obj.energy_plots(rho.get(), t.get(), n.get(), a.get()))
        all_en.grid(row=2, column=1)

        # Statistical Quantities
        # tkinter variables not converted to Python variables

        rdf = ttk.Button(text="RDF",
                         command=lambda: obj.rdf(rho.get(), t.get(), n.get(), a.get()))
        rdf.grid(row=1, column=3)

        msd = ttk.Button(text="MSD",
                         command=lambda: obj.msd(rho.get(), t.get(), n.get(), a.get()))
        msd.grid(row=1, column=4)

        vaf = ttk.Button(text="VAF",
                         command=lambda: obj.vaf(rho.get(), t.get(), n.get(), a.get()))
        vaf.grid(row=1, column=5)

        dif_plot = ttk.Button(text="D vs a",
                              command=lambda: obj.diffusion_plot(rho.get(), t.get(), n.get(), a.get()))
        dif_plot.grid(row=2, column=3)

        pc = ttk.Button(text="Pc",
                        command=lambda: obj.pc(rho.get(), t.get(), n.get(), a.get()))
        pc.grid(row=2, column=4)
        # Allows multiple figures to be stacked and then plotted
        # use tk.Button since ttk has no easy way for bg/fg manipulation
        plot_button = tk.Button(text="PLOT", bg="blue",
                                command=lambda: plt.show())
        plot_button.grid(row=1, column=7, padx=5)

        clear_figure = tk.Button(text="Clear Fig", bg="red",
                                 command=lambda: plt.clf())
        clear_figure.grid(row=1, column=8, padx=5)

        # No-Data plots
        no_data_label = ttk.Label(text="Theoretical", font=LARGE_FONT)
        no_data_label.grid(row=4, column=0, columnspan=2)
        # todo: this is a hack. find better sol'n
        empty_label = ttk.Label(text=" ")
        empty_label.grid(row=3, column=0)
        self.grid_rowconfigure(3, minsize=1000)

        # Buttons
        potential = ttk.Button(text="Potential",
                               command=lambda: obj.potential(n.get(), a.get()))
        potential.grid(row=5, column=0)

        force = ttk.Button(text="Force",
                           command=lambda: obj.force(n.get(), a.get()))
        force.grid(row=5, column=1)

        rdf2 = ttk.Button(text="RDF2",
                          command=lambda: obj.rdf_2(n.get(), a.get()))
        rdf2.grid(row=6, column=0)

        avg_q_label = ttk.Label(text="Average Quantities", font=LARGE_FONT)
        avg_q_label.grid(row=5, column=3, columnspan=3)
        three_d_label = ttk.Label(text="3D", font=LARGE_FONT)
        three_d_label.grid(row=9, column=0, columnspan=2)


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
