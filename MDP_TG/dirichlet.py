# -*- coding: utf-8 -*-

import numpy as np


class  dirichlet_dist(object):
    #----dirichlet distribution----
    def __init__(self, alpha, b):
        self.alpha = alpha
        self.b = b

    def expect(self):
        sum_alpha = sum(self.alpha)
        mean =  [a/sum_alpha for a in self.alpha]
        return mean

    def sample(self, N=1):
        s = np.random.dirichlet(self.alpha, N)
        return s




class prod_dirichlet(object):
    #----product of two dirichlet----
    def __init__(self, d1, d2):
        self.d_1 = d1
        self.d_2 = d2

    def est_mean_sigma(self, N=100):
        s1 = self.d_1.sample(N)
        s2 = self.d_2.sample(N)
        N_1 = len(self.d_1.b)
        N_2 = len(self.d_2.b)
        new_b = dict()
        for n1 in range(N_1):
            for n2 in range(N_2):
                b1 = self.d_1.[n1]
                b2 = self.d_2.[n2]
                b = (b1,b2)
                new_b[b] = []
        for i in range(N):
            for j in range(N):
                p1 = s1[i]
                p2 = s2[j]
                for n1 in range(N_1):
                    for n2 in range(N_2):
                        b1 = self.d_1.[n1]
                        b2 = self.d_2.[n2]
                        b = (b1,b2)
                        new_b[b].append(p1[n1]*p2[n2])
        mean_b = dict()
        sigma_b = dict()
        for key,value in new_b.iteritems():
            mean_p = sum(value)/(N*N)
            mean_b[key] = mean_p
            dif_v = [v-mean_p, for v in value]
            neg_dif_v = []
            for v in dif_v:
                if v>0:
                    neg_dif_v.append(0)
                else:
                    neg_dif_v.append(v)            
            sigma_b[key] = sum(neg_dif_v)/(N*N)
        return mean_b, sigma_b
                        
        
        
