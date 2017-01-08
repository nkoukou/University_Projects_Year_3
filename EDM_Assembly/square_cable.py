'''
Models the electric potential of a long square cable and analyses the numerical
solution.
'''

import numpy as np
import matplotlib.pylab as plt
import base_class as bc
reload(bc)

class SquareCable(bc.PotentialGrid):
    '''
    Represents an infinitely long square coaxial cable.
    '''
    def __init__(self, bar_width, tube_width, scale, potential):
        '''
        Input parameters:
          bar_width - width of square bar (dimensionless - integer)
          tube_width - width of square tube (dimensionless - integer)
          scale - number of grid points per unit length (integer)
          potential - fixed potential of cable (dimensionless - float)
        '''
        # Determine unit length and unit potential
        self.dim = {'x': bc.sqc_x, 'u': bc.sqc_u}

        self.scale = int(scale)
        self.bar_width = float(bar_width)
        self.tube_width = float(tube_width)
        self.u = float(potential)

        self.bar_size = int( (bar_width*scale+1)//2*2 )
        self.tube_size = int( (tube_width*scale+1)//2*2 )
        self.grid, self.fix = bc.gen_sc_grid(self.bar_size,
                                             self.tube_size, self.u)

    def set_scale(self, scale, silent=False):
        self.scale = int(scale)
        self.bar_size = int( (self.bar_width*scale+1)//2*2 )
        self.tube_size = int( (self.tube_width*scale+1)//2*2 )
        self.gen_grid()
        if not silent: print self

    def set_bar_width(self, bar_width, silent=False):
        self.bar_width = float(bar_width)
        self.bar_size = int( (bar_width*self.scale+1)//2*2 )
        self.gen_grid()
        if not silent: print self

    def set_tube_width(self, tube_width, silent=False):
        self.tube_width = float(tube_width)
        self.tube_size = int( (tube_width*self.scale+1)//2*2 )
        self.gen_grid()
        if not silent: print self

    def set_potential(self, potential, silent=False):
        self.u = float(potential)
        self.gen_grid()
        if not silent: print self

    def __repr__(self):
        return ("SquareCable( bar_width={0}, tube_width={1}, scale={2}, "\
                +"potential={3} )").format(self.bar_width, self.tube_width,
                                        self.scale, self.u)

    def gen_grid(self):
        '''
        Generates 2D system potential grid.
        '''
        self.grid, self.fix = bc.gen_sc_grid(self.bar_size,
                                             self.tube_size, self.u)

    def analyse_width_ratio(self, ratios=-1, w=-1, accuracy=0.05):
        '''
        Plots cross sections of potential for given ratios of bar to tube
        widths, relaxation parameter w and accuracy of convergence.
        '''
        if ratios==-1:
            ratios = [1.1, 1.5, 3., 6., 9., 12.]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_title(r'Potential cross section at different width ratios',
        #             fontsize=55)
        ax.set_xlabel(r'$system\ size\ (%s)$'% (self.dim['x'][1]), fontsize=40)
        ax.set_ylabel(r'$potential\ \phi\ (%s)$'% (self.dim['u'][1]),
                      fontsize=40)
        ax.tick_params(axis='both', labelsize=30)

        xvalues = int(ratios[-1]*self.bar_width*self.scale + 2)
        edge = xvalues*self.dim['x'][0]/self.scale/2.
        xaxis = np.linspace(-edge, edge, xvalues)
        for ratio in ratios:
            self.set_tube_width(ratio*self.bar_width, silent=True)
            self.converge_grid(w, accuracy)
            mid = self.grid.shape[0]/2
            xsection = np.pad(self.grid[mid]*self.dim['u'][0],
                              (xvalues/2 - mid), 'constant')
            ax.plot(xaxis, xsection, label=r'$%.1f$'% ratio)
        leg = ax.legend(prop={'size':40})
        leg.draggable()








