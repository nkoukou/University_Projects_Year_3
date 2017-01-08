'''
Defines the base class of an electric potential grid.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
from numba import jit

# Global dimensions (used for plots)
sqc_x = (2., 'cm')  # unit length for SquareCable
sqc_u = (10., 'V')  # unit potential for SquareCable
edm_x = (10., 'mm') # unit length for Edm
edm_u = (2., 'kV')  # unit potential for Edm

# Plot parameters
font = {'family' : 'normal',
        'weight' : 'normal'}
mpl.rc('font', **font)
mpl.rcParams['lines.linewidth'] = 5.

# Functions compiled just-in-time
@jit
def gen_sc_grid(b, t, u):
    '''
    Generates SquareCable grid.
    '''
    grid =  np.full((b,b), u)
    fix = np.ones((b,b))
    grid = np.pad(grid, ((t-b)/2,), 'constant', constant_values=(0,))
    fix = np.pad(fix, ((t-b)/2,), 'constant', constant_values=(0,))
    grid = np.pad(grid, 1, 'constant', constant_values=(0,))
    fix = np.pad(fix, 1, 'constant', constant_values=(1,))

    return grid, fix

@jit
def gen_edm_grid(tube_dist, scale=1):
    '''
    Generates Edm grid.
    '''
    small_plate = np.full(2,1, dtype='float64')
    big_plate = np.full(20,4, dtype='float64')
    gap = np.zeros(1, dtype='float64')

    row_one = np.concatenate((small_plate, gap, big_plate, gap, small_plate))
    row_two = np.zeros(row_one.size)
    row_three = -row_one

    grid = np.vstack((row_one, row_two, row_three))
    grid = np.pad(grid, tube_dist, 'constant', constant_values=(0,))

    fix = np.where(grid==0, 0, 1)

    if scale != 1:
        scale = np.ones((scale, scale))
        grid = np.kron(grid, scale)
        fix = np.kron(fix, scale)

    grid = np.pad(grid, 1, 'constant', constant_values=(0,))
    fix = np.pad(fix, 1, 'constant', constant_values=(1,))

    return grid, fix

@jit
def update(grid, fix, scale, w=-1):
    '''
    Updates SquareCable or Edm grid.
    Relaxation parameter w (0 < w < 2) affects the speed of convergence.
      - w = 'j': solves with Jacobi method
      - w = -1: solves with estimated optimal w
    '''
    if w=='j' or w=='J':
        new_grid=np.copy(grid)
        for index, fixed in np.ndenumerate(fix):
            if fixed: continue
            new_grid[index] = 0.25*( grid[index[0]-1, index[1]] +
                                     grid[index[0]+1, index[1]] +
                                     grid[index[0], index[1]-1] +
                                     grid[index[0], index[1]+1]  )
        return new_grid
    if w==-1:
        coef = float(grid.shape[1])/grid.shape[0]
        const = 2.0 if coef==1. else 5.5
        w = 2./(1+const/(coef*scale))


    for index, fixed in np.ndenumerate(fix):
        if fixed: continue
        grid[index] = ((1-w) * grid[index] + 0.25 * w *
                      ( grid[index[0]-1, index[1]] +
                        grid[index[0]+1, index[1]] +
                        grid[index[0], index[1]-1] +
                        grid[index[0], index[1]+1]  ))
    return grid

# Base class
class PotentialGrid(object):
    def update_grid(self, w=-1):
        '''
        Updates grid once.
        '''
        self.grid = update(self.grid, self.fix, self.scale, w)

    def converge_grid(self, w=-1, accuracy=0.05):
        '''
        Updates grid until convergence.
        '''
        temporal_spread = 1.
        spatial_spread = 0.
        updates = 0

        while temporal_spread > accuracy*spatial_spread:
            horizontal_spread = np.absolute(np.diff(self.grid, axis=-1)).max()
            vertical_spread = np.absolute(np.diff(self.grid, axis=0)).max()
            spatial_spread = max(horizontal_spread, vertical_spread)

            old_grid = np.copy(self.grid)
            self.update_grid(w)
            temporal_spread = np.linalg.norm( (self.grid - old_grid) )

            updates += 1
            if updates%1000==0:
                print '\nspatial spread = ', spatial_spread
                print 'temporal spread = ', temporal_spread
                print 'updates = ', updates
        return temporal_spread, spatial_spread, updates

    def plot_grid(self, title=None):
        '''
        Plots grid's potential field. Parameter title sets the title of the
        plot.
        '''
        if self.grid.shape[0] == self.grid.shape[1]:
            colour, shrink, aspect = 'YlOrRd', 1, (1, 10)
        else:
            colour, shrink, aspect = 'RdYlBu', 0.5, (1.2, 8)

        grid = self.dim['u'][0]*self.grid
        xedge = (grid.shape[1]-2.)*self.dim['x'][0]/self.scale/2.
        yedge = (grid.shape[0]-2.)*self.dim['x'][0]/self.scale/2.

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if title=='intro':
            ax.set_title(r'EDM experiment plate assembly', fontsize=45)
        elif title=='results':
            ax.set_title(r'Electric Potential Field', fontsize=45)
        axx = ax.imshow(grid, extent= [-xedge, xedge, -yedge, yedge],
                aspect=aspect[0], interpolation='None',
                cmap=plt.cm.get_cmap(colour))
        ax.set_xlabel(r'$system\ size\ ({0})$'.format(self.dim['x'][1]),
                      fontsize=45)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.tick_params(axis='both', labelsize=40)

        cbar = fig.colorbar(axx, shrink=shrink, aspect=aspect[1])
        cbar.ax.tick_params(labelsize=40)
        cbar.set_label(r'$Potential\ \phi\ ({0})$'.format(self.dim['u'][1]),
                       size=50)

    def analyse_scale(self, w=-1, datapoints=20, accuracy=0.05, plot=True):
        '''
        Plots number of updates against scale for given relaxation parameter w,
        number of datapoints and accuracy of convergence. If plot=False,
        returns computed updates and scales.

        Plots also maximum spatial spread of potential against scale.
        '''
        scales = np.linspace(10, 10*datapoints, datapoints)
        mesh, updates = [], []
        for s in scales:
            print s
            self.set_scale(s, silent=True)
            data = self.converge_grid(w, accuracy)
            updates.append(data[2])
            mesh.append(data[1]*self.dim['u'][0])
        if not plot: return scales, updates

        if w=='j':
            xaxis = scales*scales
            lab= r'$scale^2\ \left(\frac{1}{(%g%s)^2}\right)$'% (
                   self.dim['x'][0], self.dim['x'][1])
        else:
            xaxis = scales
            lab= r'$scale\ \left(\frac{1}{%g%s}\right)$'% (self.dim['x'][0],
                                                           self.dim['x'][1])
        slope = updates[-1]/xaxis[-1]
        fit = slope*xaxis

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(r'Number of updates against Scale', fontsize=45)
        ax.plot(xaxis, updates, label=r'Numerical data')
        ax.plot(xaxis, fit, label=r'Linear fit ($slope=%.2f$)'% (slope))
        ax.set_xlabel(lab, fontsize=35)
        ax.set_ylabel(r'$temporal\ updates$', fontsize=35)
        ax.tick_params(axis='both', labelsize=25)
        ax.legend(loc='upper left', prop={'size':40})

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(r'Spatial spread against Scale', fontsize=45)
        ax.plot(scales, mesh)
        ax.set_xlabel(r'$scale\ \left(\frac{1}{%g%s}\right)$'%
                      (self.dim['x'][0], self.dim['x'][1]), fontsize=40)
        ax.set_ylabel(r'$spatial\ spread\ (%s)$'% (self.dim['u'][1]),
                      fontsize=40)
        ax.tick_params(axis='both', labelsize=25)

    def analyse_spread(self, w=-1, datapoints=10):
        '''
        Plots spatial spread of potential against accuracy of convergence for
        given relaxation parameter w and number of datapoints.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_title(r'Spatial spread against Accuracy of convergence',
        #             fontsize=75)
        ax.set_xlabel(r'$fraction\ of\ spatial\ spread$', fontsize=40)
        ax.invert_xaxis()
        ax.set_ylabel(r'$spatial\ spread\ (%s)$'% (self.dim['u'][1]),
                      fontsize=40)
        ax.tick_params(axis='both', labelsize=30)

        accuracies = np.logspace(-1,-10,datapoints)
        for scale in np.linspace(10,10*datapoints,datapoints):
            self.set_scale(scale, silent=True)
            spreads = []
            for acc in accuracies:
                t,s,u = self.converge_grid(w, acc)
                spreads.append(s*self.dim['u'][0])
            ax.plot(accuracies, spreads, label='Scale={0}'.format(scale))

        return accuracies, spreads

    def analyse_omega(self, guess, scales=-1, datapoints=20,
                      accuracy=0.05, plot=True):
        '''
        Plots number of updates against relaxation parameter for given initial
        guess, system scales, number of datapoints and accuracy of convergence.
        If plot=False, returns computed updates and relaxation parameters.
        '''
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(r'Optimal omega search at different scales',
                         fontsize=55)
            ax.set_xlabel('$relaxation\ parameter\ \omega$', fontsize=37)
            ax.set_ylabel('$temporal\ updates$', fontsize=37)
            ax.tick_params(axis='both', labelsize=30)

        ws = np.pad(np.array([guess]), datapoints/2, 'linear_ramp',
                    end_values=(guess-0.05, 1.99))
        if scales==-1: scales = [self.scale]
        for scale in scales:
            updates = []
            for w in ws:
                self.set_scale(scale, silent=True)
                data = self.converge_grid(w, accuracy)
                updates.append(data[-1])
            if plot: ax.plot(ws, updates, label=r'Scale ${0}$'.format(scale))
            else: return ws, updates
        if plot: ax.legend(loc='upper center', prop={'size':40})

    def plot_omega_vs_scale(self, const=2., datapoints=20):
        '''
        Plots relaxation parameter against scale along with approximate fit for
        given number of datapoints.
        The fitting is approximated by the user with the constant const which
        appears in the formula: 2(1+const/(coef*scale)), where coef is the
        ratio of x and y dimensions of the system.
        '''
        coef = float(self.grid.shape[1]-2)/(self.grid.shape[0]-2)
        scales = np.linspace(10, 50, datapoints)
        fit = 2./(1+const/(coef*scales))
        ws = []
        for scale in scales:
            self.set_scale(scale, silent=True)
            guess = 2./(1+const/(coef*self.scale))
            w, update = self.analyse_omega(guess, plot=False)
            w = w[update.index(min(update))]
            ws.append(w)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(r'Relaxation parameter against scale', fontsize=55)
        ax.set_xlabel(r'$scale\ \left(\frac{1}{%g%s}\right)$'%
                      (self.dim['x'][0], self.dim['x'][1]), fontsize=37)
        ax.set_ylabel('$relaxation\ parameter\ \omega$', fontsize=37)
        ax.tick_params(axis='both', labelsize=30)
        ax.plot(scales, ws, label=r'Numerical data')
        ax.plot(scales, fit, label=r'Approximate fit ($C=%.1f$)'% (const))
        ax.legend(loc='upper left', prop={'size':40})

        return scales, ws






