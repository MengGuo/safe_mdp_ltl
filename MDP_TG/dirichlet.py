# -*- coding: utf-8 -*-

import numpy as np


class  dirichlet_dist(object):
    #----dirichlet distribution----
    def __init__(self, alpha, b):
        self.alpha = alpha
        self.b = b

    def expect(self):
        sum_alpha = sum(self.alpha)
        mean =  [a*1.0/sum_alpha for a in self.alpha]
        return mean

    def sample(self, N=1):
        s = np.random.dirichlet(self.alpha, N)
        return s


def compute_sigma(n, alpha, N=50):
    if len(alpha) <= 1:
        mean_sigma = 0
    else:
        S = np.random.beta(alpha[n], sum(alpha)-alpha[n], N)
        mean = alpha[n]*1.0/sum(alpha)
        sum_sigma = 0
        for s in S:
            if s-mean<0:
                sum_sigma += (s-mean)
        mean_sigma = sum_sigma*1.0/N
    return mean_sigma

    
def est_prod_mean_sigma(dirichlet, f_x, f_u, N=10):
    s_u_p, s_l_p = dirichlet[:]
    t_x_k = s_u_p[(f_x,f_u)]
    b_1 = t_x_k.keys()
    alpha_1 = [t_x_k[b] for b in b_1]
    d_1 = dirichlet_dist(alpha_1, b_1)
    mean_1 = d_1.expect()
    N_1 = len(b_1)
    mean_b = dict()
    sigma_b = dict()
    for n1 in range(N_1):
        sigma_1 = compute_sigma(n1, alpha_1)
        t_x = b_1[n1]        
        l_x_k = s_l_p[t_x]        
        b_2 = l_x_k.keys()
        alpha_2 = [l_x_k[b] for b in b_2]
        d_2 = dirichlet_dist(alpha_2, b_2)
        mean_2 = d_2.expect()
        N_2 = len(b_2)        
        #--------------------
        for n2 in range(N_2):
            l_x = b_2[n2]
            b = (t_x, l_x)
            sigma_2 = compute_sigma(n2, alpha_2)
            mean_b[b] = mean_1[n1]*mean_2[n2]
            sigma_b[b] = -sigma_1*sigma_2
    return mean_b, sigma_b



def est_mean_sigma(alpha, b, N=50):
    d = dirichlet_dist(alpha, b)
    S = d.sample(N)
    mean_b = dict()
    sigma_b = dict()
    for k, b_k in enumerate(b):
        mean_b[b_k] = sum([s[k] for s in S])/N        
        dif_b = [s[k]-mean_b[b_k] for s in S]
        neg_dif_b = []
        for b in dif_b:
            if b>0:
                neg_dif_b.append(0)
            else:
                neg_dif_b.append(b)
        sigma_b[b_k] = sum(neg_dif_b)/(N)
    return mean_b, sigma_b    


