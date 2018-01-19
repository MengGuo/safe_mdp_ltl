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


        
def est_mean_sigma(dirichlet, f_x, f_u, N=100):
    s_u_p, s_l_p = dirichlet[:]
    t_x_k = s_u_p[(f_x,f_u)]
    b_1 = t_x_k.keys()
    alpha_1 = [t_x_k[b] for b in b_1]
    d_1 = dirichlet_dist(alpha_1, b_1)
    s1 = d_1.sample(N)
    N_1 = len(b_1)
    out_come = dict()
    for t_x in b_1:
        l_x_k = s_l_p[t_x]
        b_2 = l_x_k.keys()
        alpha_2 = [l_x_k[b] for b in b_2]
        d_2 = dirichlet_dist(alpha_2, b_2)
        s2 = d_2.sample(N)
        N_2 = len(b_2)
        #--------------------        
        for n1 in range(N_1):
            for n2 in range(N_2):
                b1 = d_1[n1]
                b2 = d_2[n2]
                b = (b1, b2)
                outcome[b] = []
        for i in range(N):
            for j in range(N):
                p1 = s1[i]
                p2 = s2[j]
                for n1 in range(N_1):
                    for n2 in range(N_2):
                        b1 = d_1[n1]
                        b2 = d_2[n2]
                        b = (b1,b2)
                        if b not in out_come:
                            out_come[b] = [p1[n1]*p2[n2], ]
                        else:
                            out_come[b].append(p1[n1]*p2[n2])
    mean_b = dict()
    sigma_b = dict()
    for key,value in out_come.iteritems():
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
                        
        
        
