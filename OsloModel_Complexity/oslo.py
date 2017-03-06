'''
Implements the Oslo algorithm.
'''
import numpy as np
from numba import jitclass, int8, int32, float32

@jitclass([
('n', int32),
('p', float32),
('counts', int32),
('h', int32),
('slope', int8[:]),
('thresh', int8[:]),
('heights', int32[:]),
('avalanches', int32[:]),
('drops', int32[:]),
])
class System(object):
    """
    Represents a lattice which behaves according to the Oslo model.
    """
    def __init__(self, size, probability, counts):
        """
        The system size is the number of piles of the system, while the
        probability determines the slope threshold values of each pile.
        """
        self.n = int(size)
        if self.n == 1:
            raise ValueError("Algorithm not implemented for systems of size 1")
        self.p = float(probability)
        self.counts = int(counts)

        self.h = 0
        self.slope = np.zeros(self.n, dtype = np.int8)
        self.thresh = np.zeros(self.n, dtype = np.int8)

        self.heights = np.zeros(counts, dtype=np.int32)
        self.avalanches = np.zeros(counts, dtype=np.int32)
        self.drops = np.zeros(counts, dtype=np.int32)

    def initialise(self):
        '''
        Initialises the system at the empty state.
        '''
        self.h = 0
        self.slope = np.zeros(self.n, dtype = np.int8)
        for i in xrange(self.n):
            self.thresh[i] = 1 if np.random.random() < self.p else 2

    def drive(self):
        '''
        Drives the system by adding a grain on the first pile.
        '''
        self.slope[0] +=1
        self.h +=1

    def relax(self, j):
        '''
        Relaxes all piles of the system.
        '''
        while np.any(self.slope>self.thresh):
            for i in xrange(0, self.n):
                if self.slope[i] > self.thresh[i]:
                    self.avalanches[j] +=1
                    if i==0:
                        self.slope[i] -=2
                        self.slope[i+1] +=1
                        self.h -=1
                    elif i==self.n-1:
                        self.drops[j] +=1
                        self.slope[i] -=1
                        self.slope[i-1] +=1
                        self.thresh[i] = 1 if np.random.random()<self.p else 2
                    else:
                        self.slope[i] -=2
                        self.slope[i+1] +=1
                        self.slope[i-1] +=1
                    self.thresh[i] = 1 if np.random.random()<self.p else 2
        self.heights[j] += self.h

    def iterate(self):
        '''
        Iterates the algorithm.
        '''
        self.initialise()
        for j in xrange(self.counts):
            if j%int(1e5)==0: print(j)
            self.drive()
            self.relax(j)

def save(iters=1.5e7, min_size=3, max_size=11):
    '''
    Saves .npy files with raw data for heights, avalanches and grain drops.
    Runs algorithm for system sizes up to 2^max_size iterating 10^order_count
    times.
    '''
    string = lambda n: str(n).zfill(4)
    sizes = [2**x for x in range(min_size, max_size)]
    iters = int(iters)

    for s in sizes:
        print 'Size: ', s
        lat = System(s, 0.5, iters)
        lat.iterate()
        np.save('s'+string(s)+'_h', lat.heights)
        np.save('s'+string(s)+'_s', lat.avalanches)
        np.save('s'+string(s)+'_d', lat.drops)
