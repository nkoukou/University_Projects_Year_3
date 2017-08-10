'''
This module contains analysis for a growing network.
'''
import numpy as np
import scipy.optimize as sco
import scipy.stats as scs
import matplotlib as mat
import matplotlib.pylab as plt
import matplotlib.cm as cm
from graph import Graph
import logbin as lb

font = {'size'   : 20}
mat.rc('font', **font)

class Analysis(object):
    def __init__(self, N, m, iters=100, method=0, mode='temp'):
        '''
        Instatiates an object which contains methods for analysing growing
        networks.
        Parameter iters indicates how many networks will be used to improve
        statistics.
        Parameter method indicates how the network grows:
          - 0: preferential attachment to vertex proportional to its degree
          - 1xx: random attachment to vertex, after random walk of length xx
        Parameter mode indicates how the data are loaded:
          - temp: runs the networks just-in-time and directly loads the data
          - save: runs the networks just-in-time and saves the data as well
          - load: loads already saved data
        '''
        self.N = int(N)
        self.m = int(m)

        self.iters = iters
        if method in np.append(np.arange(100, 110), 0):
            self.method = method
        else:
            raise ValueError('Not implemented method')
        if mode in ['temp', 'save', 'load']:
            self.mode = mode
        else:
            raise ValueError('Not implemented mode')

        self.k = None
        self.p = None

        self.degrees = None
        self.kmax = None

        self.xb = None
        self.degb = None
        self.degb_err = None
        self.pb = None

    def setmethod(self, method):
        '''
        Setter for the method of attaching edges from a new vertex of the graph.
        '''
        if method in np.append(np.arange(100, 200), 0):
            self.method = method
        else:
            raise ValueError('Not implemented method')

    def setiters(self, iters):
        '''
        Setter for number of iterations to create new graphs.
        '''
        self.iters = int(iters)

    def setm(self, m):
        '''
        Setter for number of attachments m of new vertex.
        '''
        self.m = int(m)
        self.degrees = None
        self.kmax = None

    def setN(self, N):
        '''
        Setter for number N of final vertices in the network.
        '''
        self.N = int(N)
        self.degrees = None
        self.kmax = None

    def setkp(self):
        '''
        Setter for the range of degree k and the theoretical distribution
        across that range.
        '''
        if self.degrees is None:
            self.run_graph()
        self.k = np.arange(self.m, self.degrees.shape[1]).astype(np.float)
        if self.method==0:
            self.p = 2.*self.m*(self.m+1)/(self.k*(self.k+1)*(self.k+2))
        elif self.method==100:
            self.p = 1./(1+self.m)*np.power(self.m/(1.+self.m), self.k-self.m)

    def name_file(self, i):
        '''
        Creates the file name that is going to be loaded or saved based on
        given iteration number i.
        '''
        tail = '_N{0}_m{1}_i{2}_{3}'.format(int(np.log10(self.N)),
          str(self.m).zfill(2), str(i).zfill(3), str(self.method).zfill(3))
        if self.mode=='load':
            tail +='.npy'
        dstr = 'runs/deg'+tail
        return dstr

    def save_all(self):
        '''
        Saves all relevant data for analysis.
        '''
        if self.mode!='save':
            raise ValueError('Object not in save mode')
        self.setm(3)
        self.setiters(1000)

        for method in [0, 100]:
            self.setmethod(method)
            for N in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:
                self.setN(N)
                if self.N==int(1e6):
                    for m in [1, 3, 5, 10]:
                        self.setm(m)
                        self.run_graph()
                    self.setm(3)
                else:
                    self.run_graph()

    def run_graph(self):
        '''
        Loads data from networks.
        '''
        if self.mode=='save':
            g = Graph(self.N, self.m, self.method)
            for i in xrange(self.iters):
                g.initialize()
                g.grow()
                deg = g.deg()
                dstr = self.name_file(i)
                np.save(dstr, deg)
        else:
            self.kmax = np.zeros(self.iters, dtype=np.int64)
            self.degrees = np.zeros((self.iters,1), dtype=np.int64)
            if self.mode=='temp':
                g = Graph(self.N, self.m, self.method)
            for i in xrange(self.iters):
                if self.mode=='temp':
                    g.initialize()
                    g.grow()
                    deg = g.deg()
                elif self.mode=='load':
                    dstr = self.name_file(i)
                    deg = np.load(dstr)

                self.kmax[i] = deg.size - 1
                if deg.size > self.degrees.shape[1]:
                    padding = deg.size - self.degrees.shape[1]
                    self.degrees = np.pad(self.degrees, ((0,0),(0,padding)),
                                          'constant', constant_values=0)
                else:
                    padding = self.degrees.shape[1] - deg.size
                    deg = np.pad(deg, (0, padding), 'constant',
                                 constant_values=0)
                self.degrees[i] = deg
        self.setkp()

    def calc_stats(self, stat='chisq', notail=850):
        '''
        Analyses statistically how close the observed distribution is to the
        theoretical one.
        Parameter stat determines the statistical test used:
            - chisq ~ Pearson chi squared (plots for all tails if notail=-1)
            - redchisq ~ Reduced chi squared
            - ks ~ Kolmogorov Smirnov one-sample test
            - ks2 ~ Kolmogorov Smirnov two-sample test
        '''
        if self.degrees is None:
            self.run_graph()
        y = self.degrees[:,self.m:].sum(0)

        if stat=='chisq':
            p = self.N * self.iters * self.p
            if notail==-1:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                chis = []
                for notail in xrange(1,p.size):
                    yv, pv = y[:notail], p[:notail]
                    chi2 = scs.chisquare(yv, pv, ddof=0)
                    chis.append(chi2[1])
                ax.plot(np.arange(1, p.size)/1000., chis, 'g.', markersize=1)
                ax.set_xlabel('Degree $k/10^3$')
                ax.set_ylabel('$p$-value')
                ax.axhline(y=0.1, c='r', ls='--')
            else:
                y, p = y[:notail], p[:notail]
                chi2 = scs.chisquare(y, p, ddof=0)
                return chi2 # np.sum((np.absolute(y-p)-0.5)**2/p)
                            # - Yule's correction
        elif stat=='redchisq':
            p = self.N * self.iters * self.p
            sy = self.N * self.iters * self.degrees[:,self.m:].std(0)
            y, p, sy = y[:notail], p[:notail], sy[:notail]
            c = (y - p)/sy
            chi2 = np.sum(c*c)
            return chi2 / (notail - 1)
        elif stat=='ks':
            y = y/np.sum(y).astype(np.float)
            ks = np.random.choice(self.k, size=int(1e6), p=y)
            ks = np.sort(ks)
            indic = np.arange(1, len(ks)+1)/float(len(ks))
            cdf = np.cumsum(self.p)

            it = np.bincount(ks.astype(np.int))[self.m:]
            ix = np.cumsum(it)[:-1]
            ds = np.absolute(indic[ix-1] - cdf[:int(ks.max()-self.m)])

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.k, cdf, 'k-')
            ax.set_xlabel('Degree $k$')
            ax.set_ylabel('Theoretical $CDF$ - Observed $ECDF$')
            ax.plot(ks, indic, 'r-',
                    label='$ECDF\ d_{KS} = %f$'%(np.around(ds.max(), 5)))
            ax.legend(loc='lower right')
            return ds.max()
        elif stat=='ks2':
            y = y/np.sum(y).astype(np.float)
            ks = np.random.choice(self.k, size=int(1e6), p=y)
            ps = np.random.choice(self.k, size=int(1e6), p=self.p/np.sum(self.p))

            ks2 = scs.ks_2samp(ks[:notail], ps[:notail])
            return ks2

    def process_deg(self):
        '''
        Log bins data and creates theoretical distributions over the log-binned
        range.
        '''
        xbs = np.zeros(1)
        degbs = np.zeros((self.iters, 1))
        for i in xrange(self.iters):
            deg = self.degrees[i]
            if self.N > deg.size:
                deg = np.pad(deg, (0, self.N-deg.size), 'constant',
                             constant_values=0)
            nodes = np.repeat(np.arange(deg.size), deg).astype('int')
            xb, degb = lb.log_bin(nodes, bin_start=self.m, first_bin_width=1.,
                          a=1.25, datatype='int', drop_zeros=True)
            if xb.size > xbs.size:
                padding = xb.size - xbs.size
                xbs = xb
                degbs = np.pad(degbs, ((0,0),(0,padding)),
                                      'constant', constant_values=0)
            else:
                padding = xbs.size - xb.size
                degb = np.pad(degb, (0,padding), 'constant', constant_values=0)
            degbs[i] = degb
        self.xb = xbs
        self.degb = degbs.mean(0)
        self.degb_err = degbs.std(0)
        if self.method==0:
            self.pb = 2.*self.m*(self.m+1)/(xbs*(xbs+1)*(xbs+2))
        elif self.method==100:
            self.pb = 1./(1+self.m) * np.power(self.m/(1.+self.m), xbs-self.m)

    def plot_coefs(self, Ns=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7], ms=[1,3,5,10]):
        '''
        Plots coefficient of k_max vs N against m.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ms = np.array(ms)
        coefs = np.zeros(len(ms), dtype='float')
        coefs_err = np.zeros(len(ms), dtype='float')
        for i in xrange(len(ms)):
            print 'm =', ms[i], 'loading..'
            self.setm(ms[i])
            kavg = np.zeros(len(Ns), dtype='float')
            kerr = np.zeros(len(Ns), dtype='float')
            for j in xrange(len(Ns)):
                self.setN(Ns[j])
                self.run_graph()
                kavg[j] = self.kmax.mean()
                kerr[j] = self.kmax.std()
            f = lambda x, m, b: m*x+b
            if self.method==0:
                z, cov = sco.curve_fit(f, np.log(Ns), np.log(kavg),
                                       sigma=kerr/kavg, absolute_sigma=True)
                coefs[i] = np.exp(z[1])
                coefs_err[i] = coefs[i]*np.sqrt(cov[1][1])
            elif self.method==100:
                z, cov = sco.curve_fit(f, np.log(Ns), kavg,
                                       sigma=kerr, absolute_sigma=True)
                y = 1./np.log(1+1./ms)
                coefs[i] = z[0]
                coefs_err[i] = np.sqrt(cov[0][0])
        if self.method==0:
            y = np.sqrt(ms*(ms+1.))
        if self.method==100:
            y = 1./np.log(1+1./ms)

        ax.errorbar(ms, coefs, yerr=coefs_err, fmt='bo', label='Observed')
        ax.plot(ms, y, 'r--', label='Theory')
        ax.set_xlabel('Edges added per step $m$')
        ax.set_ylabel('$k_{max}\ coefficient\ C$')
        ax.legend(loc='upper left')


    def plot_ms(self, ms=[1,3,5,10]):
        '''
        Plots degree distribution for different m.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        c = cm.rainbow(np.linspace(0, 1, len(ms)))
        for i in xrange(len(ms)):
            print 'm =', ms[i], 'loading..'
            self.setm(ms[i])
            self.run_graph()
            self.process_deg()
            xb, degb, degb_err = self.xb, self.degb, self.degb_err
            if len(ms)==1:
                clr = 'r'
                deg = self.degrees.sum(0) / float(self.degrees.sum())
                ax.loglog(self.k, deg[self.m:], 'b.')
            else:
                clr = c[i]
            errax = ax.errorbar(xb, degb, yerr=degb_err, color=clr, fmt='-',
                                lw=1.5, capsize=2, errorevery=1,
                                label='$m = {0}$'.format(ms[i]))
            errax[-1][0].set_linestyle('-')
            errax[-1][0].set_linewidth(0.8)
            ax.loglog(self.k, self.p, 'k--', lw=1.2)
        ax.set_xlabel('Degree $k$')
        ax.set_ylabel('Degree distribution $p(k)$')
        ax.legend(loc='lower left')

    def plot_Ns(self, Ns=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7]):
        '''
        Plots distribution for different N, k_max against N and collapses data.
        '''
        Ns = np.array(Ns, dtype='float')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        figkn = plt.figure()
        axkn = figkn.add_subplot(111)
        figcol = plt.figure()
        axcol = figcol.add_subplot(111)
        c = cm.rainbow(np.linspace(0, 1, len(Ns)))
        kavg = np.zeros(len(Ns), dtype='float')
        kerr = np.zeros(len(Ns), dtype='float')
        for i in xrange(len(Ns)):
            print 'N =', int(np.log10(Ns[i])), 'loading..'
            self.setN(Ns[i])
            self.run_graph()
            kavg[i] = self.kmax.mean()
            kerr[i] = self.kmax.std()

            self.process_deg()
            xb, degb, degb_err = self.xb, self.degb, self.degb_err
            errax = ax.errorbar(xb, degb, yerr=degb_err, color=c[i], fmt='-',
                                lw=1.5, capsize=2, errorevery=1,
                                label='$10^{0}$'.format(int(np.log10(Ns[i]))))
            errax[-1][0].set_linestyle('-')
            errax[-1][0].set_linewidth(0.8)
            ax.loglog(self.k, self.p, 'k--')
            axcol.loglog(xb/kavg[i], degb/self.pb, color=c[i],
                         label='$10^{0}$'.format(int(np.log10(Ns[i]))))
        ax.set_xlabel('Degree $k$')
        ax.set_ylabel('Degree distribution $p(k)$')
        ax.legend(loc='upper right', ncol=1, prop={'size':20})
        axcol.set_xlabel('Scaled degree $k/k_{max}$')
        axcol.set_ylabel('Scaled distribution $p(k)/P(k)$')
        axcol.axhline(y=1, c='k', ls='--')
        axcol.legend(loc='lower left', ncol=2)

        f = lambda x, m, b: m*x+b
        if self.method==0:
            z, cov = sco.curve_fit(f, np.log(Ns), np.log(kavg),
                                       sigma=kerr/kavg, absolute_sigma=True)
            y = np.exp(z[1])*Ns**z[0]
        elif self.method==100:
            z, cov = sco.curve_fit(f, np.log(Ns), kavg,
                                       sigma=kerr, absolute_sigma=True)
            y = z[1] + z[0]*np.log(Ns)
        axkn.errorbar(Ns, kavg, yerr=kerr, fmt='bo')
        axkn.plot(Ns, y, 'r--')
        axkn.set_xlabel('System size $N$')
        axkn.set_ylabel('Maximum degree $k_{max}$')
        return z[0], z[1], cov

    def plot_ls(self, els=[0, 1, 2, 5, 10]):
        '''
        Performs analysis on different lengths of random walker, producing
        relevant plots.
        '''
        els = np.array(els)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        figsl = plt.figure()
        axsl = figsl.add_subplot(111)
        c = cm.rainbow(np.linspace(0, 1, len(els)))

        zs = np.zeros(len(els), dtype='float')
        for i in xrange(len(els)):
            print 'el =', els[i], 'loading..'
            self.setmethod(100+els[i])
            self.run_graph()
            self.process_deg()
            xb, degb = self.xb, self.degb
            ax.loglog(xb, degb, '-', color=c[i], lw=1.5,
                        label='$\ell = {0}$'.format(els[i]))

            z = np.polyfit(np.log(xb), np.log(degb), 1)
            zs[i] = z[0]

        for method in [0, 100]:
            self.setmethod(method)
            self.setkp()
            ax.loglog(self.k, self.p, 'k--')
        ax.set_xlabel('Degree $k$')
        ax.set_ylabel('Degree distribution $p(k)$')
        ax.legend(loc='upper right', ncol=2, prop={'size':20})

        if self.m==1:
            axsl.plot(els[els%2==0], zs[els%2==0], 'x-.', c[0], lw=1.5,
                      markersize=9)
            axsl.plot(els[els%2==1], zs[els%2==1], 'x-.', c[-1], lw=1.5,
                      markersize=9)
        else:
            axsl.plot(els, zs, 'rx-.', lw=1.5, markersize=9)
        axsl.axhline(y=-3, c='k', ls='--', lw=1.)
        axsl.set_xlabel('Walk length $\ell$')
        axsl.set_ylabel('Degree distribution slope $\gamma$')










