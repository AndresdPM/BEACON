"""
This is the OPTICS core for BEACON
"""

import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.spatial.distance as H
import sys

class tree_node(object):
    def __init__(self, points, start, end, parentNode):
        self.points = points
        self.start = start
        self.end = end
        self.parentNode = parentNode
        self.children = []
        self.splitpoint = -1

    def __str__(self):
        return "start: %d, end %d, split: %d" % (self.start, self.end, self.splitpoint)

        
    def assignSplitPoint(self,splitpoint):
        self.splitpoint = splitpoint

    def addChild(self, child):
        self.children.append(child)

def optics(x, k, distMethod = 'euclidean'):
    if len(x.shape)>1:
        m,n = x.shape
    else:
        m = x.shape[0]
        n == 1

    try:
        D = H.squareform(H.pdist(x, distMethod))
        distOK = True
    except Exception, ex:
        print ex
        print "Squareform or pdist error"
        distOK = False

    CD = np.zeros(m)
    RD = np.ones(m)*1E10

    for i in xrange(m):
        tempInd = D[i].argsort()
        tempD = D[i][tempInd]
        CD[i] = tempD[k]#**2

    order = []
    seeds = np.arange(m, dtype = np.int)

    ind = 0
    while len(seeds) != 1:
        ob = seeds[ind]
        seedInd = np.where(seeds != ob)
        seeds = seeds[seedInd]

        order.append(ob)
        tempX = np.ones(len(seeds))*CD[ob]
        tempD = D[ob][seeds]#[seeds]

        temp = np.column_stack((tempX, tempD))
        mm = np.max(temp, axis = 1)
        ii = np.where(RD[seeds]>mm)[0]
        RD[seeds[ii]] = mm[ii]
        ind = np.argmin(RD[seeds])


    order.append(seeds[0])
    RD[0] = 0
    return RD, CD, order

def euclid(i, x):
    y = np.zeros_like(x)
    y += 1
    y *= i
    if len(x) != len(y):
        raise ValueError, "Vectors must be same length!"

    d = (x-y)**2
    return np.sqrt(np.sum(d, axis = 1))

def local_maxima(index,Reach_Plot,Reach_Points,nghsize):

    for i in range(1,nghsize+1):
        if index + i < len(Reach_Plot):
            if (Reach_Plot[index] < Reach_Plot[index+i]):
                return 0
            
        if index - i >= 0:
            if (Reach_Plot[index] < Reach_Plot[index-i]):
                return 0
    
    return 1

def get_leaves(node, arr):
    if node is not None:
        if node.splitpoint == -1:
            arr.append(node)
        for n in node.children:
            get_leaves(n,arr)
    return arr

def find_local_maxima(Reach_Plot, Reach_Points, nghsize):
    
    localMaximaPoints = {}
    
    for i in range(1,len(Reach_Points)-1):
        if Reach_Plot[i] > Reach_Plot[i-1] and Reach_Plot[i] >= Reach_Plot[i+1] and local_maxima(i,Reach_Plot,Reach_Points,nghsize) == 1:
            localMaximaPoints[i] = Reach_Plot[i]
    
    return sorted(localMaximaPoints, key=localMaximaPoints.__getitem__ , reverse=True)

def get_array(node,num, arr):
    if node is not None:
        if len(arr) <= num:
            arr.append([])
        try:
            arr[num].append(node)
        except:
            arr[num] = []
            arr[num].append(node)
        for n in node.children:
            get_array(n,num+1,arr)
        return arr
    else:
        return arr
 
def cluster_tree(node, parentNode, localMaximaPoints, Reach_Plot, Reach_Points, min_cluster_size, maximaRatio = 0.75):
    if len(localMaximaPoints) == 0:
        return
    
    s = localMaximaPoints[0]
    node.assignSplitPoint(s)
    localMaximaPoints = localMaximaPoints[1:]

    Node1 = tree_node(Reach_Points[node.start:s],node.start,s, node)
    Node2 = tree_node(Reach_Points[s+1:node.end],s+1, node.end, node)
    LocalMax1 = []
    LocalMax2 = []

    for i in localMaximaPoints:
        if i < s:
            LocalMax1.append(i)
        if i > s:
            LocalMax2.append(i)
    
    Nodelist = []
    Nodelist.append((Node1,LocalMax1))
    Nodelist.append((Node2,LocalMax2))
    
    significantMin = .003

    if Reach_Plot[s] < significantMin:
        node.assignSplitPoint(-1)
        cluster_tree(node,parentNode, localMaximaPoints, Reach_Plot, Reach_Points, min_cluster_size, maximaRatio)
        return

    checkRatio = .8
    checkValue1 = int(np.round(checkRatio*len(Node1.points)))
    checkValue2 = int(np.round(checkRatio*len(Node2.points)))
    if checkValue2 == 0:
        checkValue2 = 1
    avgReachValue1 = float(np.average(Reach_Plot[(Node1.end - checkValue1):Node1.end]))
    avgReachValue2 = float(np.average(Reach_Plot[Node2.start:(Node2.start + checkValue2)]))


    rejectionRatio = maximaRatio-0.05

    if float(avgReachValue1 / float(Reach_Plot[s])) > maximaRatio or float(avgReachValue2 / float(Reach_Plot[s])) > maximaRatio:

        if float(avgReachValue1 / float(Reach_Plot[s])) < rejectionRatio:
            Nodelist.remove((Node2, LocalMax2))
        if float(avgReachValue2 / float(Reach_Plot[s])) < rejectionRatio:
            Nodelist.remove((Node1, LocalMax1))
        if float(avgReachValue1 / float(Reach_Plot[s])) >= rejectionRatio and float(avgReachValue2 / float(Reach_Plot[s])) >= rejectionRatio:
            node.assignSplitPoint(-1)
            cluster_tree(node,parentNode, localMaximaPoints, Reach_Plot, Reach_Points, min_cluster_size, maximaRatio)
            return
 
    #remove clusters that are not big enought
    if len(Node1.points) < min_cluster_size:
        try:
            Nodelist.remove((Node1, LocalMax1))
        except Exception:
            sys.exc_clear()
    if len(Node2.points) < min_cluster_size:
        try:
            Nodelist.remove((Node2, LocalMax2))
        except Exception:
            sys.exc_clear()
    if len(Nodelist) == 0:
        node.assignSplitPoint(-1)
        return
    

    similaritythreshold = 0.4
    bypassNode = 0
    if parentNode != None:
        sumRP = np.average(Reach_Plot[node.start:node.end])
        sumParent = np.average(Reach_Plot[parentNode.start:parentNode.end])
        if float(float(node.end-node.start) / float(parentNode.end-parentNode.start)) > similaritythreshold:
            parentNode.children.remove(node)
            bypassNode = 1
        
    for nl in Nodelist:
        if bypassNode == 1:
            parentNode.addChild(nl[0])
            cluster_tree(nl[0], parentNode, nl[1], Reach_Plot, Reach_Points, min_cluster_size, maximaRatio)
        else:
            node.addChild(nl[0])
            cluster_tree(nl[0], node, nl[1], Reach_Plot, Reach_Points, min_cluster_size, maximaRatio)

def auto_cluster(Reach_Plot, Reach_Points, min_cluster_size_ratio, maxima_ratio = 0.75):

    min_neighborhood_size = 2
    min_maxima_ratio = 0.001
    
    min_cluster_size = int(min_cluster_size_ratio * len(Reach_Points))
    
    nghsize = int(min_maxima_ratio*len(Reach_Points))

    if nghsize < min_neighborhood_size:
        nghsize = min_neighborhood_size
    
    localMaximaPoints = find_local_maxima(Reach_Plot, Reach_Points, nghsize)
    
    rootNode = tree_node(Reach_Points, 0, len(Reach_Points), None)
    cluster_tree(rootNode, None, localMaximaPoints, Reach_Plot, Reach_Points, min_cluster_size, maxima_ratio)

    return rootNode

"""
Andres del Pino Molina
"""
