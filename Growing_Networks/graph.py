'''
This module contains the class of a growing network.
'''
from numba import jitclass, int32
import numpy as np

@jitclass([
('N', int32),
('m', int32),
('random', int32),
('edg_num', int32),
('edges', int32[:]),
])
class Graph(object):
    '''
    Represents a growing graph.
    '''
    def __init__(self, N, m, method):
        self.N = int(N)
        self.m = int(m)
        self.random = int(method)

        self.edg_num = 2*self.m*self.N - self.m*(self.m+1)
        self.edges = np.zeros(self.edg_num, dtype=np.int32)

    def initialize(self):
        '''
        Initializes a complete graph with (self.m) vertices.
        '''
        self.edges = np.zeros(self.edg_num, dtype=np.int32)
        i = 0
        for m in xrange(self.m):
            for n in xrange(m+1, self.m):
                self.edges[2*i] = m+1
                self.edges[2*i+1] = n+1
                i +=1

    def pick_stub(self, n):
        '''
        Picks vertex to connect to new vertex according to given method.
        '''
        ind = 2*self.m*n - self.m*(self.m+1)
        if self.m==1 and ind==0:
            return np.array([1], dtype=np.int32)
        stubs = np.zeros(self.m, dtype=np.int32)

        if self.random>100:
            length = self.random-100
            for i in xrange(self.m):
                stub = int(np.random.uniform(1, n+1))
                lg = 0
                while lg < length:
                    nns = np.where(self.edges[:ind]==stub)[0]
                    index = int(np.random.uniform(0, nns.size))
                    vertex = nns[index]
                    if self.edges[:ind][vertex-1] < self.edges[vertex]:
                        stub = self.edges[vertex-1]
                    else:
                        stub = self.edges[vertex+1]
                    lg +=1
                stubs[i] = stub

        elif self.random==100:
            for i in xrange(self.m):
                stub = int(np.random.uniform(1, n+1))
                while np.any(stubs==stub):
                    stub = int(np.random.uniform(1, n+1))
                stubs[i] = stub
        else:
            for i in xrange(self.m):
                index = int(np.random.uniform(0, ind))
                while np.any(stubs==self.edges[index]):
                    index = int(np.random.uniform(0, ind))
                stubs[i] = self.edges[index]
        return stubs

    def add_node(self, n):
        '''
        Adds a new node to the graph.
        '''
        ind = 2*self.m*n - self.m*(self.m+1)
        stubs = self.pick_stub(n)
        for i in xrange(self.m):
            self.edges[ind+2*i] = stubs[i]
            self.edges[ind+2*i+1] = n+1

    def grow(self):
        '''
        Adds nodes until graph has self.N number of vertices.
        '''
        for n in xrange(self.m, self.N):
            self.add_node(n)

    def deg(self):
        '''
        Returns degree distribution from edge list.
        '''
        nodes = np.bincount(self.edges)[1:].astype(np.int32)
        return np.bincount(nodes)



