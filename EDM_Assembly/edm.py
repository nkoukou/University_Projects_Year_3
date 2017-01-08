'''
Models the EDM experiment electric plate assembly.
'''

import numpy as np
import matplotlib.pylab as plt
import base_class as bc
reload(bc)

class Edm(bc.PotentialGrid):
    '''
    Represents the EDM experiment electric field plate setup.
    '''
    def __init__(self, tube_dist, scale):
        '''
        Input parameters:
          tube_dist - distance to the grounded border (dimensionless - integer)
          scale - number of grid points per unit length (1/[10mm] - integer)
        '''
        self.dim = {'x': bc.edm_x, 'u': bc.edm_u}

        self.tube_dist = tube_dist
        self.scale = scale

        self.grid, self.fix = bc.gen_edm_grid(tube_dist, scale)

    def set_tube_dist(self, tube_dist, silent=False):
        self.tube_dist = int(tube_dist)
        self.gen_grid()
        if not silent: print self

    def set_scale(self, scale, silent=False):
        self.scale = int(scale)
        self.gen_grid()
        if not silent: print self

    def __repr__(self):
        return "Edm( tube_dist={0}, scale={1} )".format(self.tube_dist,
                self.scale)

    def gen_grid(self):
        '''
        Generates 2D grid for EDM experiment.
        '''
        self.grid, self.fix = bc.gen_edm_grid(self.tube_dist, self.scale)

    def calc_efield(self, plot=True):
        '''
        Calculates and plots electric field along particle path.
        '''
        h = self.dim['x'][0]/self.scale
        mid = self.grid.shape[0]/2
        path = self.grid[mid-2:mid+3,:]
        efield = []
        for i in range(1, path.shape[1]-1):
            efx = (+ path[2,i+1] - path[2,i-1])*self.dim['u'][0]/(2*h)
            efy = ((- path[0,i] + 8*path[1,i] - 8*path[3,i] + path[4,i]
                   )*self.dim['u'][0])/(12*h)
            efield.append(np.sqrt(efx*efx + efy*efy))
        path = path[:,1:-1]

        edge = len(path[2])*self.dim['x'][0]/self.scale/2.
        xaxis = np.linspace(-edge, edge, len(path[2]))
        if not plot: return xaxis, efield

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_title(r'Cross section of electric field along particle path',
        #             fontsize=55)
        ax.set_xlabel(r'$system\ size\ (%s)$'% (self.dim['x'][1]), fontsize=37)
        ax.set_ylabel(r'$electric\ field\ E\ (MVm^{-1})$', fontsize=37)
        ax.tick_params(axis='both', labelsize=27)
        ax.plot(xaxis, efield)

    def analyse_homogeneity(self, tube_dist=1, plot=True):
        '''
        Calculates homogeneity of electric field along particle path for given
        tube distance from the system.
        '''
        self.set_tube_dist(tube_dist, silent=True)
        self.converge_grid(w=-1, accuracy=5e-11)
        xaxis, efield = self.calc_efield(plot=False)
        small_homogen, big_homogen, big_ind = [], [], (tube_dist+3)*self.scale

        small_plate = efield[tube_dist*self.scale+1:(tube_dist+2)*self.scale]
        big_plate = efield[big_ind+1:-big_ind]
        small_ref = small_plate[len(small_plate)/2]
        big_ref = big_plate[len(big_plate)/2]
        small_check = np.absolute(small_plate - small_ref)/small_ref
        big_check = np.absolute(big_plate - big_ref)/big_ref

        levels = np.logspace(0,-5,1000)
        for level in levels:
            small_pass = small_check[small_check<level]
            big_pass = big_check[big_check<level]
            small_homogen.append(len(small_pass))
            big_homogen.append(len(big_pass))

        small_homogen = np.array(small_homogen)*self.dim['x'][0]/self.scale
        big_homogen = np.array(big_homogen)*self.dim['x'][0]/self.scale

        if not plot: return levels, small_homogen, big_homogen

        fig = plt.figure()
        ax = fig.add_subplot(111)
        axb, axs = fig.add_subplot(211), fig.add_subplot(212)

        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top='off', bottom='off',
                       left='off', right='off')
        ax.set_title(r'Electric field homogeneity between plates', fontsize=55)
        ax.set_xlabel(r'$level\ of\ homogeneity$', fontsize=37, labelpad=5)
        ax.set_ylabel(r'$path\ length\ (%s)$'% (self.dim['x'][1]), fontsize=37,
                      labelpad=20)

        axb.plot(levels, big_homogen, label=r'Big plates')
        axb.legend(loc='lower left', prop={'size':40})
        axb.tick_params(axis='both', labelsize=27)
        axb.invert_xaxis()
        axs.plot(levels, small_homogen, color='g', label=r'Small plates')
        axs.legend(loc='lower left', prop={'size':40})
        axs.tick_params(axis='both', labelsize=27)
        axs.invert_xaxis()

    def dustify(self, plate='central', pos='centre', behaviour='conductor'):
        '''
        Adds a dust particle of size 100um at given plate, on side facing the
        particle path. Can be positioned (pos) at the middle or at the central
        edge of plate (left edge for central plate).

        It can either be an insulator or a conductor corresponding to a fixed
        potential of 0 or the value at given plate respectively.
        '''
        row = 1+(self.tube_dist+1)*self.scale
        if plate=='left':
            if pos=='centre': col = (self.tube_dist+1)*self.scale
            else: col = (self.tube_dist+2)*self.scale
        elif plate=='central':
            if pos=='centre': col = (self.tube_dist+13)*self.scale
            else: col = 1+(self.tube_dist+3)*self.scale
        elif plate=='right':
            if pos=='centre': col = (self.tube_dist+25)*self.scale
            else: col = 1+(self.tube_dist+24)*self.scale

        if behaviour=='insulator':
            u = 0.
        elif behaviour=='conductor':
            u = self.grid[row-1, col]

        if self.scale==100:
            self.grid[row, col] = u
            self.fix[row, col] = 1.
        elif self.scale==200:
            self.grid[row:row+2, col:col+2] = u
            self.fix[row:row+2, col:col+2] = 1.
        else:
            raise ValueError('Method dustify requires scale = 100 or 200')

    def plot_middle_dust(self):
        '''
        Plots electic field for the cases of a pure system and of an
        insulating and a conducting dust particle sitting at the middle of the
        central plate.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(r'Cross section of electric field along particle path',
                     fontsize=55)
        ax.set_xlabel(r'$system\ size\ (%s)$'% (self.dim['x'][1]), fontsize=37)
        ax.set_ylabel(r'$electric\ field\ E\ (MVm^{-1})$', fontsize=37)
        ax.tick_params(axis='both', labelsize=27)

        self.converge_grid(accuracy=1e-7)
        xaxis, efield = self.calc_efield(plot=False)
        ax.plot(xaxis, efield, label=r'No dust')
        self.dustify(behaviour='conductor')
        self.converge_grid(accuracy=1e-7)
        xaxis, efield = self.calc_efield(plot=False)
        ax.plot(xaxis, efield, label=r'Conducting dust')
        self.dustify(behaviour='insulator')
        self.converge_grid(accuracy=1e-7)
        xaxis, efield = self.calc_efield(plot=False)
        ax.plot(xaxis, efield, label=r'Insulating dust')
        leg = ax.legend(prop={'size':40})
        leg.draggable()


















