'''
This module contains one function which can be used to reproduce the figures in
the report. In order to reduce time, much fewer datapoints and smaller systems
are considered compared to the report. The user can change this.

There are three more modules:
    - oslo.py (contains oslo algorithm)
    - log_bin_CN_2016.py (contains the given log-binning function)
    - analysis.py (contains a class that analyses the oslo model)
All modules have sufficient documentation for the user to explore them.
'''
import analysis as an
reload(an)

def plot_figure(number, counts=int(1e6), sizes=an.all_sizes[:6]):
    '''
    Plots figure which has given number in the report.
    Parameter counts determines iterations of the model.
    Parameter sizes determines the sizes to be considered.
    '''
    lat = an.Analysis(counts, sizes)
    print 'Results will differ from report results due to'
    print 'low number of system sizes and datapoints used'
    if number==1:
        print '\nThis figure is referenced in its caption.'
    if number==2:
        lat.test_btw(size=64)
    if number==3:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        print lat.havg, lat.tc, lat.sizes
        lat.plot_heights()
    if number==4:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        lat.plot_crossover()
    if number==5:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        lat.collapse_heights()
    if number==6:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        lat.scale_havg()
    if number==7:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        lat.scale_hstd()
    if number==8:
        lat.get_data(1, 'auto')
        lat.calc_havg(); lat.calc_tc(); lat.calc_hstd()
        lat.prob_height(False)
        lat.prob_height(True)
    if number==9:
        lat.get_data(2, 'auto')
        lat.prob_aval(int(1e6))
    if number==10:
        lat.get_data(2, 'auto')
        lat.scale_aval(False, False)
        lat.scale_aval(True, False)
    if number==11:
        lat.get_data(2, 'auto')
        lat.moments()
    if number==12:
        lat.get_data(3, 'auto')
        lat.scale_aval(False, True)
        lat.scale_aval(True, True)