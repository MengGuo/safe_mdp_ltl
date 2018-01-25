from math import sqrt, atan2

from random import random

def distance(x1, x2):
    return sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2)

    
class sensor(object):

    def __init__(self, real_map, radius, decay, slope):
        self.real_map = real_map
        self.radius = radius
        self.decay = decay
        self.slope = slope

        def measure_by_sense(prod_mdp, current_state):
            # measurement from the sensor
            # s_p={(x,u,x'):k,}
            # l_p={(x,l):k,}
            #
            s_p = dict()
            l_p = dict()            
            f_x = current_state[0]
            # for height difference
            # --------------------
            for x in prod_mdp.nodes():
                dist_xy = distance(f_x, x)
                if (dist_xy < self.raidus):
                    # height
                    real_height = self.real_map.node[x]['height']
                    meas_height = real_height *(1 + randam()*self.decay*dist_xy)
                    prod_mdp.graph['mdp'].node[x]['height'] = meas_height
                    prod_mdp.node[x]['height'] = meas_height
            # --------------------
            f_height = prod_mdp.graph['mdp'].node[x]['height'] 
            for t_x in prod_mdp.graph['mdp'].successors_iter(f_x):
                dist_xy = distance(f_x, t_x)
                t_height = prod_mdp.graph['mdp'].node[t_x]['height']
                dif_height = abs(t_height-f_height)
                if (dist_xy >0) and (slope[0]<=atan2(dif_height, dist_xy)<=slope[1]):
                    prop = prod_mdp.graph['mdp'].edge[f_x][t_x]
                    for u in prop.keys():
                        s_p[(f_x, u, t_x)] = -10
                        


            return s_p, l_p
