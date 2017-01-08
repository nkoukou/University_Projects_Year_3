'''
Produces the figures included in the report.
'''

import numpy as np
import matplotlib.pylab as plt

import square_cable as sq
import edm as e
reload(e)
reload(sq)

def plot_figure(number):
    '''
    Plots figure given the number it has in the report.
    '''
    print '\nFontsizes correspond to report size.\n'
    print 'Maximising the figure, zooming in'
    print 'or converting to log scales may be helpful.'

    if number==1:
        sys = e.Edm(1,100)
        sys.plot_grid(title='intro')
    if number==2:
        sys = sq.SquareCable(1,3,100,1)
        sys.plot_grid(title=None)
    if number==3:
        sys = sq.SquareCable(1,3,20,1)
        sys.analyse_spread(w=-1, datapoints=10)
    if number==4:
        sys = sq.SquareCable(1,3,100,1)
        sys.converge_grid(w=-1, accuracy=0.05)
        sys.plot_grid(title=None)
    if number==5:
        sys = sq.SquareCable(1,3,10,1)
        sys.analyse_scale(w='j', datapoints=15, accuracy=0.05, plot=True)
    if number==6:
        sys = sq.SquareCable(1,3,10,1)
        sys.analyse_scale(w=-1, datapoints=20, accuracy=0.05, plot=True)
    if number==7:
        sys = sq.SquareCable(1,3,20,1)
        sys.analyse_width_ratio(ratios=-1, w=-1, accuracy=0.05)
    if number==8:
        sys = e.Edm(1,100)
        sys.converge_grid(w=-1, accuracy=0.05)
        sys.plot_grid(title=None)
    if number==9:
        sys = e.Edm(1,100)
        sys.converge_grid(w=-1, accuracy=0.05)
        sys.calc_efield(plot=True)
    if number==10:
        sys = e.Edm(1,100)
        sys.analyse_homogeneity(tube_dist=1, plot=True)
    if number==11:
        sys = e.Edm(1,100)
        sys.plot_middle_dust()
    if number==12:
        sys = sq.SquareCable(1,3,20,1)
        sys.analyse_omega(guess=1.8, scales=[60,80,100], datapoints=20,
                          accuracy=0.05, plot=True)
    if number==13:
        sys = sq.SquareCable(1,3,20,1)
        sys.plot_omega_vs_scale(const=2., datapoints=20)
    if number==14:
        sys = e.Edm(1,100)
        sys.plot_omega_vs_scale(const=5.5, datapoints=20)

###############################################################################
##Uncomment next line and vary argument by figure number in the report (1-14)##
#plot_figure(number=2)


