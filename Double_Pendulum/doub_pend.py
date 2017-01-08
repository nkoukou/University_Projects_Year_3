'''
Solves numerically and analyses the stability of the solution for the
double pendulum using scaled, dimensionless variables.
'''

import numpy as np
import matplotlib.pylab as plt

class DoubleBob(object):
    '''
    Represents the moving bobs in a double pendulum. Object B hangs from
    object A.
    All attributes are private and the user is advised to set them only via the
    setter methods provided, but underscores have been omitted for clarity.
    '''
    def __init__(self, t_f, R, G, y0, h):
        '''
        Parameters:
          t_f - scaled duration of experiment
          R - mass ratio, A over B
          G - scaled damping coefficient
          y0 - initial angular amplitudes and velocities of objects A and B.
               In particular, y0 = (theta, phi, theta-dot, phi-dot), where
               theta and phi are object A and B 's amplitudes respectively
          h - scaled time step
        All parameters are numbers, except y0 which is a list of four numbers.
        '''
        self.tf = float(t_f)
        self.R = float(R)
        self.G = float(G)
        self.y0 = np.array(y0, dtype=float)
        self.h = float(h)

        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0

    def set_tf(self, tf, silent=False):
        self.tf = float(tf)
        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0
        if not silent: print self

    def set_R(self, R, silent=False):
        self.R = float(R)
        if not silent: print self

    def set_G(self, G, silent=False):
        self.G = float(G)
        if not silent: print self

    def set_y0(self, y0, silent=False):
        self.y0 = np.array(y0, dtype=float)
        self.y[:,0] = self.y0
        if not silent: print self

    def set_h(self, h, silent=False):
        self.h = float(h)
        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0
        if not silent: print self

    def __repr__(self):
        """
        Provides information on DoubleBob object.
        """
        return "DoubleBob(t_f = {0}, R = {1}, G = {2}, y0 = {3}, h = {4})"\
               .format(self.tf, self.R, self.G, self.y0, self.h)

    def runge_kutta(self):
        '''
        Solves system using the Runge-Kutta 4 method.
        '''
        h, R, G = self.h, self.R, self.G
        f = np.array([ [0, 0, 1, 0],
                       [0, 0, 0, 1],
                       [-(R+1), R, -G, 0],
                       [(R+1), -(R+1), G*(1-1/R), -G/R] ])
        for i in range(1, self.t.size):
            k1 = f.dot(self.y[:,i-1])
            y2 = self.y[:,i-1] + 0.5*h*k1
            k2 = f.dot(y2)
            y3 = self.y[:,i-1] + 0.5*h*k2
            k3 = f.dot(y3)
            y4 = self.y[:,i-1] + 1.*h*k3
            k4 = f.dot(y4)
            self.y[:,i] = self.y[:,i-1] + h/6. * (k1 + 2*k2 + 2*k3 + k4)

    def energy_ratio(self, plot_energy=True):
        '''
        Returns and plots energy over initial energy.
        Potential energy is set to 0 at the lowest level the pendulum reaches.
        '''
        y, R = self.y, self.R
        T = 0.5*y[2]*y[2] + 0.5*R*y[2]*y[2] + 0.5*R*y[3]*y[3] + \
            R*y[2]*y[3]*np.cos(y[0] - y[1])
        V = -(1+R)*np.cos(y[0]) - R*np.cos(y[1]) + (2*R+1)
        E = T + V

        if plot_energy:
            plt.plot(self.t, E/E[0])
            plt.xlabel('Time')
            plt.ylabel('Fractional energy')
        return E/E[0]


    def stab_analysis(self, dps=3):
        '''
        Analyses stability of pendulum solutions and determines time step at
        which the energy blows up.
        Parameter dps determines the number of accurate decimal places in the
        solution.
        '''
        # Set reasonable values for h and tf to speed up the algorithm
        self.set_h(0.01, silent=True)
        self.set_tf(500, silent=True)

        # Set up the bisection method to find critical h
        upper_bound = 1.1
        h_prev = 3.
        diff = h_prev-self.h

        # Run the algorithm
        self.runge_kutta()
        ratio = self.energy_ratio(plot_energy=False)[-1]
        stable = ratio < upper_bound

        while stable and (diff>10**(-dps)):
            self.set_h(0.5*(self.h+h_prev), silent=True)
            diff = h_prev-self.h
            self.runge_kutta()
            ratio = self.energy_ratio(plot_energy=False)[-1]
            stable = ratio < upper_bound
            while not stable:
                h_prev = self.h
                self.set_h(self.h-diff/2., silent=True)
                self.runge_kutta()
                ratio = self.energy_ratio(plot_energy=False)[-1]
                stable = ratio < upper_bound
        print self
        return self.h