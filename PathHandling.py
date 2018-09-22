import platform
import os
"""
This is legacy code that should not be used, in any case!
"""

class OSPaths:
    """
    Managing paths in Windows and Linux
    need to import os and platform
    """
    # __pathWIN = 'C:/Users/user/Desktop/MD simulation/Archives of Data/'
    # __pathLINUX = '/media/gn/ECD68B8AD68B53AC/Users/user/Desktop/MD simulation/Archives of Data/'
    # This is not a Valid Path anymore, repair disc
    # Always Launch from MD simulation dir
    # __pathLINUX = '/media/gn/4C1EE3431EE324A4/Code/C++/MD simulation/Archives of Data/'
    __pathWIN = '../Archives of Data/tests/4^3/'
    __pathLINUX = __pathWIN
    __iso = '/Isothermal~'
    __non_iso = '/Non-Isothermal'
    __density = ['Density 0.5', 'Density 0.8', 'Density 1.2', 'Density 3.6']
    __step = 'step '
    __steps_num = ['50', '2500', '3000', '5000', '10000', '12500', '15000', '20000', '5000 v2']

    def __init__(self):
        self.OS = platform.system()
        if self.OS == 'Linux':
            self.path = self.__pathLINUX
            print(self.OS)                   # self.path is now
        else:                               # publicly available
            self.path = self.__pathWIN
            print(self.OS)

    def dir(self, density, steps):
        if density not in self.__density:
            raise ValueError("Not a valid Density.")
        if steps not in self.__steps_num:
            raise ValueError("Not a valid number of steps.")

        den = density
        ste = steps
        directory = self.path + den + self.__iso + self.__step + ste
        os.chdir(directory)
        return directory

