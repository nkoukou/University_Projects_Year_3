'''
Analyses the Oslo algorithm.
'''
import numpy as np
import matplotlib as m
import matplotlib.pylab as plt
import matplotlib.ticker as mtk
import matplotlib.cm as cm
import log_bin_CN_2016 as lb
import oslo as o

font = {'size'   : 19}
m.rc('font', **font)

all_sizes = np.array([8, 16, 32, 64, 128, 256, 512, 1024, 2048])

class Analysis(object):
    """
    Includes plots, and analysis tools for a system that behaves according to
    the Oslo model..
    """
    def __init__(self, counts, sizes=all_sizes):
        '''
        The number of counts given defines the interval over which averages
        are calculated.
        '''
        self.counts = int(counts)
        self.sizes = sizes

    def set_counts(self, counts):
        self.counts = int(counts)

# Get data
    def get_data(self, data, generate='ready'):
        '''
        Gets data to be analysed. Parameter generate is set to 'ready' to
        load existing data files or 'auto' to generate data for sizes in
        sizes with iterations equal to (size^2 + self.counts) and threshold
        probability 0.5. Parameter data is 1 for heights, 2 for avalanches and
        3 for grain drops
        '''
        if data not in [1, 2, 3]:
            raise ValueError('data = 1, 2 or 3')
        for size in self.sizes:
            print size
            if generate=='ready':
                string = lambda n: str(n).zfill(4)
                if data==1:
                    y = np.load('s'+string(size)+'_h.npy')
                if data==2:
                    y = np.load('s'+string(size)+'_s.npy')
                if data==3:
                    y = np.load('s'+string(size)+'_d.npy')
            elif generate=='auto':
                lat = o.System(size, 0.5, self.counts)
                lat.iterate()
                if data==1:
                    y = lat.heights
                if data==2:
                    y = lat.avalanches
                if data==3:
                    y = lat.drops
            if size==self.sizes[0]:
                self.y = y
            else:
                self.y = np.vstack((self.y, y))
        self.x = np.arange(len(self.y[0]))

# Calculate quantities
    def moving_avg(self, data, w=25):
        '''
        Performs a moving average smoothing over the given data.
        '''
        ww= 2*w+1
        mov_avg = np.cumsum(data, dtype=float)
        mov_avg[ww:] =  (mov_avg[ww:] - mov_avg[:-ww])
        return mov_avg[ww-1:]/ww

    def calc_havg(self):
        '''
        Calculates height average after cross over time
        '''
        self.havg = np.around(self.y[:, -self.counts:].mean(1), 1)

    def calc_hstd(self):
        '''
        Calculates height standard deviationat after cross over time
        '''
        self.hstd = np.around(self.y[:, -self.counts:].std(1), 2)

    def calc_tc(self):
        '''
        Estimates cross over time
        '''
        tc = []
        for i in range(len(self.sizes)):
            diff = []
            heights = self.y[i, :2*self.sizes[i]*self.sizes[i]]
            section = self.sizes[i]/2
            splits = len(heights)/section
            for arr in np.split(heights, splits):
                diff.append(self.havg[i] - arr.mean())
            tc.append(int(np.argmin(np.array(diff)>1) + 0.5)*section)
        self.tc = np.array(tc)

    def calc_moment(self, i, k, drops):
        '''
        Calculates k-th moment of avalanches in given interval.
        '''
        if drops:
            y= self.y[i,:]
            y = y[np.nonzero(y)]
        else:
            y = self.y[i, -self.counts:]
        return 1./y.size * np.sum(np.power(y, float(k)))

    def calc_prob(self, size, counts=None, drops=True):
        '''
        Calculates pribability of heights, avalanches or grain drops.
        '''
        if counts is None: counts = self.counts
        i = np.where(self.sizes==size)[0][0]
        if drops:
            y= self.y[i,:]
            y = y[np.nonzero(y)]
        else:
            y = self.y[i, -counts:]
        return np.bincount(y)/float(y.size)

# Plots and analysis
    def test_btw(self, size, p=None):
        '''
        Tests consistency of Oslo model for given size.
        '''
        if p is None: p = [0., 0.25, 0.5, 0.75, 1.]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("$sites$")
        ax.set_ylabel("$lattice\ height\ h\ /L$")
        c = cm.rainbow(np.linspace(0, 1, len(p)))
        x = np.arange(size)
        heights = []
        for j in xrange(len(p)):
            lat = o.System(size, p[j], size*size + self.counts)
            lat.iterate()
            s = lat.slope
            h = np.zeros(len(s)); h[0] = lat.h
            heights.append(lat.h/64.)
            for i in range(1, len(s)):
                h[i] = h[i-1] - s[i-1]
            ax.plot(x, h/64., color=c[j], lw=2.,
                    label='p={0:.2f}'.format(p[j]))
        ax.legend(loc='upper right', prop={'size':17})
        return heights

    def plot_heights(self):
        '''
        Plots height against time for all system sizes.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("$number\ of\ grains\ \log_2{t}$")
        ax.set_ylabel("$height\ \log_2{h}$")
        c = cm.rainbow(np.linspace(0, 1, len(self.sizes)))
        h = np.log2(self.havg)
        t = np.log2(self.tc)
        for i in xrange(len(self.sizes)):
            x, y = np.log2(self.x), np.log2(self.y[i,:])
            ax.plot(x, y, '-', lw=0.8, color = c[i], label=self.sizes[i])
            ax.axhline(h[i], color='k', lw=0.3, linestyle='--')
            ax.axvline(t[i], color='k', lw=0.3, linestyle='--')
        leg = ax.legend(loc='upper left', ncol=2, prop={'size':15})
        for lg in leg.legendHandles:
            lg.set_linewidth(1.5)

    def plot_crossover(self):
        '''
        Plots dependence of average height and crossover time on system size.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        z = np.polyfit(self.sizes, self.havg, 1)
        ax.set_xlabel("$system\ size\ L$")
        ax.set_ylabel("$average\ height$")
        ax.loglog(self.sizes, self.havg, 'rx', markersize=20.)
        ax.loglog(self.sizes, z[0]*self.sizes+z[1], 'b-', lw=2.,
                  label='slope {0:.3f}'.format(z[0]))
        ax.legend(loc='upper left', prop={'size':15})

        fig = plt.figure()
        ss = self.sizes*self.sizes
        z = np.polyfit(ss, self.tc, 1)
        ax = fig.add_subplot(111)
        ax.set_xlabel("$system\ size\ L^2$")
        ax.set_ylabel("$crossover\ time\ t_c$")
        ax.loglog(ss, self.tc, 'rx', markersize=20.)
        ax.loglog(ss, z[0]*ss+z[1], 'b-', lw=2.,
                  label='slope {0:.3f}'.format(z[0]))
        ax.legend(loc='upper left', prop={'size':15})

    def collapse_heights(self, w=25):
        '''
        Plots the collapsed heights against time.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r"$scaled\ time\ t/L^2$")
        ax.set_ylabel(r"$scaled\ height\ h/L$")
        c = cm.rainbow(np.linspace(0, 1, len(self.sizes)))
        for i in xrange(len(self.sizes)):
            heights = self.moving_avg(self.y[i,:], w=w)
            s = float(self.sizes[i])
            x = self.x[w:int(10*self.sizes[i]**2)]/(s*s)
            y = heights[:int(10*self.sizes[i]**2-w)]/s
            ax.loglog(x, y, color=c[i], lw=1.5, label=int(s))
        ax.legend(loc='upper left', ncol=2, prop={'size':13})

        fig = plt.figure()
        ax = fig.add_subplot(111)
        tc = self.tc[-1]
        x, y = np.log2(self.x[1:tc]), np.log2(self.y[-1, 1:tc])
        z = np.polyfit(x, y, 1)
        ax.plot(x, y, 'gx')
        ax.plot(x, z[0]*x+z[1], 'k-', label='slope %g'%(z[0]))
        ax.legend(loc='upper left')
        return z[0], np.exp(z[1]), tc/(s*s), self.havg[-1]/s

    def scale_havg(self, a=None):
        '''
        Plots average height with corrections.
        '''
        if a is None: a = np.linspace(1.7328, 1.7343, 5)
        fig = plt.figure()
        fig.suptitle("Height against Size with Corrections")
        for i in range(len(a)):
            ax = fig.add_subplot(2,4,i+1)
            x = np.log(self.sizes)
            y = np.log(1-self.havg/(a[i]*self.sizes))
            z = np.polyfit(x, y, 1)
            ax.plot(x, y)
            ax.plot(x, z[0]*x+z[1])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("$system\ size\ L$")
        ax.set_ylabel("$transformed\ height\ h'$")
        xx, yy = self.sizes, 1-self.havg/(a.mean()*self.sizes)
        x = np.log2(xx)
        y = np.log2(yy)
        z = np.polyfit(x, y, 1)
        ax.plot(x, y, 'rx', markersize=15.)
        ax.plot(x, z[0]*x+z[1], 'b-')

        a0, a1 = np.around(a.mean(), 4), np.around(np.exp(z[1]), 4)
        w1 = np.around(-z[0], 4)
        return a0, a1, w1

    def scale_hstd(self):
        '''
        Plots scaling of height standard deviation.
        '''
        x, y = self.sizes, self.hstd

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("$system\ size\ L/10^3$")
        ax.set_ylabel("$standard\ deviation\ \sigma_h$")
        ax.plot(x/1e3, y, lw=2.)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        x, y = np.log10(x), np.log10(y)
        z = np.polyfit(x, y, 1)
        ax.plot(x, y, 'rx', markersize=15.)
        ax.plot(x, z[0]*x+z[1], 'b-', label='slope {0:.2}'.format(z[0]))
        ax.set_xlabel("$system\ size\ \log{L}$")
        ax.set_ylabel("$standard\ deviation\ \log{\sigma_h}$")
        ax.legend(loc='upper left', prop={'size':16})

        fig = plt.figure()
        ax = fig.add_subplot(111)
        x, y = self.sizes, self.hstd/self.sizes**z[0]
        ax.loglog(x, y, 'rx')
        ax.axhline(y=0.5838, color='k', linestyle='--')
        ax.set_xlabel("$L$")
        ax.set_ylabel("$\sigma_h/L^{0.24}$")

        return z[0], np.exp(z[1])

    def avg_slope(self):
        '''
        Plots average slope and its std.
        '''
        zavg, zstd = self.havg/self.sizes, self.hstd/self.sizes

        fig = plt.figure()
        ax = fig.add_subplot(121)
        ax.axhline(y=1.733, color='k', linestyle='--')
        ax.plot(self.sizes, zavg)
        ax = fig.add_subplot(122)
        ax.plot(self.sizes, zstd)

    def prob_height(self, collapse=True):
        '''
        Plots Gaussian height probabilities before or after collapse.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        c = cm.rainbow(np.linspace(0, 1, len(self.sizes)))
        if collapse:
            ax.set_xlabel("$scaled\ height\ z$")
            ax.set_ylabel('$probability\ \tilde{P}(z;L)$')
        else:
            ax.set_xlabel("$height\ h$")
            ax.set_ylabel('$probability\ P(h;L)$')
        for i in range(len(self.sizes)):
            prob = self.calc_prob(self.sizes[i])
            x, y = np.arange(len(prob)), prob
            #x = np.nonzero(prob)[0]
            #y = prob[x]
            if collapse:
                x, y = (x-self.havg[i])/self.hstd[i], self.hstd[i]*y
                #if self.sizes[i]==2048: x, y = 0, 0
                ax.plot(x, y, color=c[i], lw=2., label=self.sizes[i])
            else:
                ax.semilogx(x, y, color=c[i], lw=2., label=self.sizes[i])
        ax.legend(loc='upper right', ncol=2, prop={'size':13})

    def prob_aval(self, counts, base=1.2, size=all_sizes[5], drops=False):
        '''
        Plots the log-binning of avalanche probability.
        '''
        prob = self.calc_prob(size=size, counts=counts, drops=drops)
        fig = plt.figure()
        x = np.nonzero(prob)[0]
        y = prob[x]
        j = np.where(self.sizes==size)[0][0]
        if base is None: base = np.linspace(1.15, 1.25, 4)
        ax = fig.add_subplot(111)
        if drops:
            s= self.y[j,:]
            s = s[np.nonzero(s)]
            ax.set_xlabel('$drop\ size\ d$')
            ax.set_ylabel('$probability\ P(d;L)$')
        else:
            s = self.y[j, -counts:]
            ax.set_xlabel('$avalanche\ size\ s$')
            ax.set_ylabel('$probability\ P(s;L)$')
        xx, yy = lb.log_bin(s, 0., 1., base, 'integer', False)
        ax.loglog(x, y, 'x')
        ax.loglog(xx, yy, 'o', lw=2., label='a='+str(base))
        ax.legend(loc='upper right')

    def scale_aval(self, collapse=True, drops=False):
        '''
        Plots probability of avalanches or drops before or after collapse for
        all sizes.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        c = cm.rainbow(np.linspace(0, 1, len(self.sizes)))
        for i in xrange(len(self.sizes)):
            print self.sizes[i]
            if drops:
                y= self.y[i,:]
                y = y[np.nonzero(y)]
                k = [1.002, 1.249]
            else:
                y = self.y[i, -self.counts:]
                k = [1.55, 2.25]
            xx, yy = lb.log_bin(y, 0, 1., 1.2, 'integer', False)
            xx, yy = np.array(xx), np.array(yy)
            if collapse:
                y = (xx**k[0])*yy
                x = xx/(self.sizes[i]**k[1])
                ax.loglog(x, y, '-', lw=1.5, color=c[i], label=self.sizes[i])
                ax.legend(loc='lower left', ncol=2, prop={'size':13})
                if drops:
                    ax.set_xlabel('$d/L^{1.25}$')
                    ax.set_ylabel('$d^{1.01} P(d;L)$')
                else:
                    ax.set_xlabel('$s/L^{2.24}$')
                    ax.set_ylabel('$ s^{1.55} P(s;L)$')
            else:
                #x, y = np.log2(xx), np.log2(yy)
                ax.loglog(xx, yy, '-', lw=1.5, color=c[i], label=self.sizes[i])
                ax.legend(loc='upper right', ncol=2, prop={'size':13})
                if drops:
                    ax.set_xlabel('$drop\ size\ d$')
                    ax.set_ylabel('$probability\ P(d;L)$')
                else:
                    ax.set_xlabel('$avalanche\ size\ s$')
                    ax.set_ylabel('$probability\ P(s;L)$')


    def moments(self, points=5, check=5, drops=False):
        '''
        Performs moment analysis for avalanche or drop probability.
        '''
        moments = np.zeros((len(self.sizes), points))
        coefs, ycept = [], []
        for k in xrange(points):
            for i in xrange(len(self.sizes)):
                moments[i,k] = self.calc_moment(i, k+1, drops=drops)
            z = np.polyfit(np.log(self.sizes[4:]), np.log(moments[4:,k]), 1)
            coefs.append(z[0])
            ycept.append(np.exp(z[1]))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        k = check-1
        y = moments[:,k]/self.sizes**(coefs[k])
        ax.plot(self.sizes, y, 'rx', markersize=15., label='M'+str(k+1))
        ax.axhline(y=ycept[k], color='k', linestyle='--')
        ax.set_xlabel("$system\ size\ L$")
        ax.set_ylabel(r"$s^k/L^{1+k-{\tau}_s}$")
        ax.legend(loc='upper right', prop={'size':15})

        ks = np.arange(1, points+1)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        z = np.polyfit(ks, coefs, 1)
        ax.plot(ks, coefs, 'x', markersize=15.)
        ax.plot(ks, z[0]*ks+z[1], '-', lw=2., label='D = %.3g'%(z[0]))
        ax.set_xlabel("$moment\ order\ k$")
        ax.set_ylabel("$Dk+D(1-\tau_s)$")
        ax.legend(loc='upper left', prop={'size':15})

        d = z[0]
        tau = 1. - z[1]/d

        return d, tau





