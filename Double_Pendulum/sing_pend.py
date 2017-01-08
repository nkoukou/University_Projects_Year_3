'''
Solves numerically and analyses the stability of the solution for the
simple pendulumusing scaled, dimensionless variables.
'''

import numpy as np
import matplotlib.pylab as plt

class Bob(object):
    '''
    Represents the moving object in a single pendulum.
    All attributes are private and the user is advised to set them only via the
    setter methods provided, but underscores have been omitted for clarity.
    '''
    def __init__(self, t_f, G, y0, h):
        '''
        Parameters:
          t_f - scaled duration of experiment
          G - scaled damping coefficient
          y0 - initial angular amplitude and velocity, y0 = (theta, theta-dot)
          h - scaled time step
        All parameters are numbers, except y0 which is a list of two numbers.
        '''
        self.tf = float(t_f)
        self.G = float(G)
        self.y0 = np.array(y0, dtype=float)
        self.h = float(h)

        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0

    def set_tf(self, tf, silent=False):
        self.tf = tf
        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0
        if not silent: print self

    def set_G(self, G, silent=False):
        self.G = G
        if not silent: print self

    def set_y0(self, y0, silent=False):
        self.y0 = np.array(y0, dtype=float)
        self.y[:,0] = self.y0
        if not silent: print self

    def set_h(self, h, silent=False):
        self.h = h
        self.t = np.linspace(0, self.tf, int(self.tf/self.h))
        self.y = np.zeros((self.y0.size, self.t.size)); self.y[:,0] = self.y0
        if not silent: print self

    def __repr__(self):
        """
        Provides information on Bob object.
        """
        return "Bob(t_f = {0}, G = {1}, y0 = {2}, h = {3})".format(
                self.tf, self.G, self.y0, self.h)

    def explicit_euler(self):
        '''
        Solves system using the explicit Euler method.
        '''
        f = np.array([ [1, self.h], [-self.h, (1-self.h*self.G)] ])
        for i in range(1, self.t.size):
            self.y[:,i] = f.dot(self.y[:,i-1])

    def leapfrog(self):
        '''
        Solves system using the leapfrog method. To get the second angle of
        motion, an explicit Euler method with 1/10 of the original time step is
        implemented.
        '''
        h = 0.1*self.h
        tf = self.t[1]
        t = np.linspace(0, tf, int(tf/h))
        y = np.zeros((self.y0.size, t.size)); y[:,0] = self.y0
        f = np.array([ [1, h], [-h, (1-h*self.G)] ])
        for i in range(1, t.size):
            y[:,i] = f.dot(y[:,i-1])
        self.y[:,1] = y[:,-1]

        f = np.array([ [0, 2*self.h], [-2*self.h, -2*self.h*self.G] ])
        for i in range(2, self.t.size):
            self.y[:,i] = self.y[:,i-2] + f.dot(self.y[:,i-1])

    def runge_kutta(self):
        '''
        Solves system using the Runge-Kutta 4 method.
        '''
        f = np.array([ [0, 1], [-1, -self.G] ])
        for i in range(1, self.t.size):
            k1 = f.dot(self.y[:, i-1])
            y2 = self.y[:,i-1] + 0.5*self.h*k1
            k2 = f.dot(y2)
            y3 = self.y[:,i-1] + 0.5*self.h*k2
            k3 = f.dot(y3)
            y4 = self.y[:,i-1] + 1.*self.h*k3
            k4 = f.dot(y4)
            self.y[:,i] = self.y[:,i-1] + self.h/6. * (k1 + 2*k2 + 2*k3 + k4)

    def implicit_euler(self):
        '''
        Solves system using the implicit Euler method.
        '''
        u, v = self.y
        for i in range(1, self.t.size):
            v[i] = (v[i-1] - self.h*u[i-1])/(1+self.h*self.G+self.h*self.h)
            u[i] = u[i-1] + self.h*v[i]
        self.y[0], self.y[1] = u, v

    def exact_euler(self):
        '''
        Solves system using the explicit Euler method, but not the small
        angle approximation.
        '''
        u, v = self.y
        for i in range(1, self.t.size):
            v[i] = -self.h*np.sin(u[i-1]) + (1-self.h*self.G)*v[i-1]
            u[i] = u[i-1] + self.h*v[i-1]
        self.y[0], self.y[1] = u, v

    def method(self, meth='rk'):
        '''
        Chooses method to solve the system. Parameter meth can be:
          ee - explicit_euler
          lp - leapfrog
          rk - runge_kutta
          ie - implicit_euler
          ce - exact_euler
        '''
        if meth=='ee':
            self.explicit_euler()
        if meth=='lp':
            self.leapfrog()
        if meth=='rk':
            self.runge_kutta()
        if meth=='ie':
            self.implicit_euler()
        if meth=='ce':
            self.exact_euler()

    def energy_ratio(self, meth='rk', plot_energy=True):
        '''
        Returns and can plot energy over initial energy for given method.
        Potential energy is set to 0 at the lowest level the pendulum reaches.
        '''
        self.method(meth)
        y = self.y
        E = 0.5*y[1]*y[1] + (1 - np.cos(y[0]))
        if plot_energy:
            plt.plot(self.t, E/E[0])
            plt.xlabel('Time')
            plt.ylabel('Fractional energy')
        return E/E[0]

    def is_stable(self, ratio, meth, low_bounds=False):
        '''
        Returns criteria for stability.
        Parameter low_bounds checks for energy decay in an undamped case.
        Leapfrog and exact euler present big fluctuations.
        '''
        if meth=='lp': M=100.
        elif meth=='ce': M=1.e15
        else: M = 1.1
        m = 0.9

        if low_bounds:
            return ratio > m
        return ratio < M

    def stab_analysis(self, meth='rk', dps=3):
        '''
        Analyses stability of system and returns whether it is unconditionally
        stable, unstable or conditionally stable along with the threshold of
        time step for stability.
        Determines if energy blows up and at what time step, but also warns
        about unexpected energy decay in the undamped cases.
        Parameter dps is number of decimal places counted for accuracy.
        '''
        # Set reasonable values for h and tf to speed up the algorithm
        if meth=='ce':
            self.set_h(4, silent=True)
            self.set_tf(100000, silent=True)
        else:
            self.set_h(0.1, silent=True)
            self.set_tf(1000, silent=True)

        # Set up the bisection method to find critical h
        hmin = 5e-4
        if meth=='ce': hmax = 1000
        else: hmax = 10.
        h_prev = hmax
        diff = abs(h_prev-self.h)

        # Run algorithm
        ratio = self.energy_ratio(meth, plot_energy=False)[-1]
        stable = self.is_stable(ratio, meth)

        while not stable and self.h>=hmin:
            self.set_h(0.1*self.h, silent=True)
            ratio = self.energy_ratio(meth, plot_energy=False)[-1]
            stable = self.is_stable(ratio, meth)

        while stable and (self.h<=hmax) and (diff>10**(-dps)):
            self.set_h(0.5*(self.h+h_prev), silent=True)
            diff = abs(h_prev-self.h)
            ratio = self.energy_ratio(meth, plot_energy=False)[-1]
            stable = self.is_stable(ratio, meth)
            while not stable:
                h_prev = self.h
                self.set_h(self.h-diff/2., silent=True)
                ratio = self.energy_ratio(meth, plot_energy=False)[-1]
                stable = self.is_stable(ratio, meth)

        #  Check for energy decay just before blow-up
        if stable and self.G==0.:
            hi = self.h
            self.set_h(self.h-10.*hmin, silent=True)
            ratio = self.energy_ratio(meth, plot_energy=False)[-1]
            stable = self.is_stable(ratio, meth, low_bounds=True)
            if not stable:
                print "\nEnergy decays unexpectedly."
            elif stable:
                print "\nEnergy does not decay."
            self.set_h(hi, silent=True)

        # Return results of stability analysis
        if self.h<hmin:
            print "\nUnconditionally unstable. Energy blows up.\n"
            print self
            return None
        if self.h>(hmax-10.*diff):
            print "\nUnconditionally stable.\n"
            print self
            return None
        if diff<10**(-dps):
            print '\nConditionally stable under h_crit = {0}\n'.format(self.h)
            print self
            return self.h