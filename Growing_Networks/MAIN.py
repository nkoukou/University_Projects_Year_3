'''
This module reproduces the figures of the report.

The other three modules included in the project are:
    - logbin.py ~ Contains a logbinning function adapted by James Clough's code
    - graph.py ~ Contains a growing network class
    - analysis.py ~ Contains analysis methods for the growing networks

Figures produced in this module are of lower quality compared to the report,
to ensure that figures are produced fast and do not take up significant amounts
of memory.

The user is free to explore all modules, which contain sufficent documentation,
but can produce all figures with the plot_figure function below by only using
the last line of this code.
'''
from analysis import Analysis

def plot_figure(number):
    '''
    Parameter number corresponds to the number of the figure in the report.
    '''
    if number==1:
        a = Analysis(1e5, 1, 100, 0, 'temp')
        a.plot_ms(ms=[1])
        a.plot_ms(ms=[1,3,5,10])
    if number==2:
        a = Analysis(1e5, 1, 100, 0, 'temp')
        a.calc_stats(notail=-1)
    if number==3:
        a = Analysis(1e5, 1, 100, 0, 'temp')
        a.calc_stats('ks', notail=50)
    if number in [4,5,6]:
        a = Analysis(1e5, 1, 200, 0, 'temp')
        a.plot_Ns(Ns=[1e2,1e3,1e4,1e5])
        a = Analysis(1e5, 1, 100, 0, 'temp')
        a.plot_coefs(Ns=[1e2, 1e3, 1e4, 1e5], ms=[1,3,5,10])
    if number==7:
        a = Analysis(1e5, 1, 100, 100, 'temp')
        a.plot_ms(ms=[1])
        a.plot_ms(ms=[1,3,5,10])
    if number in [8,9,10]:
        a = Analysis(1e5, 1, 400, 100, 'temp')
        a.plot_Ns(Ns=[1e2,1e3,1e4,1e5])
        a = Analysis(1e5, 1, 100, 100, 'temp')
        a.plot_coefs(Ns=[1e2, 1e3, 1e4, 1e5], ms=[1,3,5,10])
    if number==11:
        a = Analysis(1e4, 2, 20, 100, 'temp')
        a.plot_ls(els=[0, 1, 2, 5, 10])
    if number==12:
        print 'Ignore solid lines on graph.'
        a = Analysis(1e4, 1, 20, 100, 'temp')
        a.plot_ls(els=[0, 1, 2, 3, 4])

#plot_figure(1)