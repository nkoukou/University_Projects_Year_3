'''
Produces useful figures for analysis of the single and double pendulum
numerical solutions and their stability.
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import sing_pend as sp
import double_pend as dp

meths = {'ee': 'Explicit Euler', 'ie': 'Implicit Euler',
         'ce': 'Exact Euler', 'lp': 'Leapfrog', 'rk': 'Runge Kutta 4'}

def plot_sing_osc(bob, methods=meths, fit=True):
    '''
    Plots oscillations of a single pendulum bob for given methods (e.g.
    ['ee', 'ie']). Parameter fit plots a cosine along with the solutions.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for meth in methods:
        bob.method(meth)
        ax.plot(bob.t, bob.y[0], lw= 5, label=meths[meth])
    if fit:
        ax.plot(bob.t, bob.y0[0]*np.cos(bob.t), lw= 5, label='Cosine fit')
    ax.tick_params(axis='both', labelsize=20)
    fig.suptitle('Comparison of numerical solutions', fontsize=45)
    plt.subplots_adjust(top=0.85)
    ax.set_xlabel(r'Time $t\ \left(\sqrt{\frac{l}{g}}\right)$', fontsize=30)
    ax.set_ylabel(r'Amplitude $\theta\ (radians)$', fontsize=30)
    ax.legend(bbox_to_anchor=(0., 1.0, 1., .10), loc=1, ncol=3,
              mode="expand", borderaxespad=0., prop={'size':30})

def plot_sing_stab(bob, ready_data):
    '''
    Plots critical time step against damping coefficient for given bob in the
    cases of conditional stability.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)

    Gsteps = np.logspace(-1, 1, 50)
    if ready_data:
        hee = np.load('ee.npy')
        hrk = np.load('rk.npy')
        hce = np.load('ce.npy')
    else:
        hee, hrk, hce = [], [], []
        for G in Gsteps:
            print 'G = ', G
            bob.set_G(G, silent=True)
            hee.append(bob.stab_analysis('ee'))
            hrk.append(bob.stab_analysis('rk'))
            hce.append(bob.stab_analysis('ce'))

    ax.plot(Gsteps, hee, lw=5., label=meths['ee'])
    ax.plot(Gsteps, hrk, lw=5., label=meths['rk'])
    ax.plot(Gsteps, hce, lw=5., label=meths['ce'])

    ax.tick_params(axis='both', labelsize=20)
    fig.suptitle('Stability dependence on damping', fontsize=45)
    plt.subplots_adjust(top=0.92)
    ax.set_xlabel(r'Damping coefficient $D\ \left(m\sqrt{gl}\right)$',
                      fontsize=30)
    ax.set_ylabel(
      r'Critical time step $h_{c}\ \left(\sqrt{\frac{l}{g}}\right)$',
      fontsize=30)
    ax.legend(loc='upper right', prop={'size':30})

def plot_doub_osc(bob):
    '''
    Plots oscillation of a double pendulum bob.
    '''
    bob.runge_kutta()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    axs, axb = fig.add_subplot(211), fig.add_subplot(212)

    bob.set_R(0.01, silent=True)
    bob.runge_kutta()
    axs.plot(bob.t, bob.y[0], color='g', label=r'$\theta$')
    axs.plot(bob.t, bob.y[1], color='b', label=r'$\phi$')
    bob.set_R(100, silent=True)
    bob.runge_kutta()
    axb.plot(bob.t, bob.y[0], color='g', label=r'$\theta$')
    axb.plot(bob.t, bob.y[1], color='b', label=r'$\phi$')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off',
                   left='off', right='off')
    ax.set_title(r'Numerical solution of motion', fontsize=45)
    plt.subplots_adjust(top=0.95)
    ax.set_xlabel(r'Time $t\ \left(\sqrt{\frac{l}{g}}\right)$', fontsize=30,
                   labelpad=5)
    ax.set_ylabel(r'Amplitudes $\theta,\ \phi\ (radians)$', fontsize=30,
                  labelpad=40)
    axs.legend(loc='upper right', prop={'size':30})
    axs.tick_params(axis='both', labelsize=20)
    axb.legend(loc='upper right', prop={'size':30})
    axb.tick_params(axis='both', labelsize=20)

def plot_doub_stab(bob, ready_data):
    '''
    Plots critical time step against mass ratio and damping coefficient for
    given double bob.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if ready_data:
        r = np.load('r.npy')
        g = np.load('g.npy')
        h = np.load('h.npy')
        h = h.reshape((len(r), len(h)/len(r)))
    else:
        Gsteps, Rsteps = np.logspace(-1, 1, 50), np.logspace(-2, 2, 50)
        hsteps = []
        for G in Gsteps:
            print 'G = ', G
            bob.set_G(G, silent=True)
            for R in Rsteps:
                print 'R = ', R
                bob.set_R(R, silent=True)
                hsteps.append(bob.stab_analysis())
        g, r = np.meshgrid(Gsteps, Rsteps)
        h = np.array(hsteps).reshape((len(Gsteps), len(hsteps)/len(Gsteps)))
    rr, gg = np.log10(r), np.log10(g)
    surf = ax.plot_surface(rr, gg, h, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig.suptitle('Stability dependence on mass ratio and damping', fontsize=45)
    plt.subplots_adjust(top=0.92)
    ax.set_xlabel(r'Mass ratio $R$',
                      fontsize=30, labelpad=25)
    ax.set_ylabel(r'Damping coefficient $D\ \left(m\sqrt{gl}\right)$',
                      fontsize=30, labelpad=25)
    ax.set_zlabel(
      r'Critical time step $h_{c}\ \left(\sqrt{\frac{l}{g}}\right)$',
                  fontsize=30, labelpad=15)
    xticks = [1e-2, 1e-1, 1e0, 1e1, 1e2]
    ax.set_xticks(np.log10(xticks))
    ax.set_xticklabels(xticks, fontsize=20)
    yticks = [1e-1, 0.316, 1e0, 3.16, 1e1]
    ax.set_yticks(np.log10(yticks))
    ax.set_yticklabels(yticks, fontsize=20)
    zmin = round(h.ravel().min(), 2)
    zmax = round(h.ravel().max(), 2)
    zticks = np.linspace(zmin, zmax, 4)
    ax.set_zticks(zticks)
    ax.set_zticklabels(zticks, fontsize=20)
    fig.colorbar(surf, shrink=0.5, aspect=5).ax.tick_params(labelsize=25)

#############
# MAIN RUNS #
#############

def plot_figure(fig, ready_data=False):
    '''
    Plots figures of the report.
    Parameter fig is the number of the figure in the report.
    Parameter ready_data, when True, skips all calculations and plot graph
    using the .npy files provided.
    '''
    if fig==1:
        bob = sp.Bob(t_f=1000, G=0., y0=[0.1, 0], h=0.001)
        plot_sing_osc(bob, methods=meths, fit=True)
    if fig==2:
        bob = sp.Bob(t_f=100, G=0., y0=[0.1, 0], h=0.1)
        plot_sing_stab(bob, ready_data)
    if fig==3:
        bob = dp.DoubleBob(t_f=100, R=1, G=0., y0=[0.1, 0, 0, 0], h=0.001)
        plot_doub_osc(bob)
    if fig==4:
        bob = dp.DoubleBob(t_f=100, R=1, G=0, y0=[0.1, 0, 0, 0], h=0.1)
        plot_doub_stab(bob, ready_data)

def plot_table(method):
    '''
    Returns data of the table of the report for given method.
    Parameter method is a key of the meths dictionary.
    '''
    print 'Method:', method
    bob = sp.Bob(100, 0, [0.1, 0], 0.1)
    bob.stab_analysis(method)
    print '-----'
    bob = sp.Bob(100, 0.2, [0.1, 0], 0.1)
    bob.stab_analysis(method)