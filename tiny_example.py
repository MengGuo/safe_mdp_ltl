from MDP_TG.mdp import Motion_MDP
from MDP_TG.dra import Dra, Product_Dra,  execution_with_sensing
from MDP_TG.vis import visualize_world_paths
from MDP_TG.sense import sensor

from bootstrap import construct_nodes, construct_edges


from math import radians

import time

t0 = time.time()

#-------- construct model -------
l = 1 #m
N = 4
label_set = set([frozenset(['o',]), frozenset(['h',]), frozenset(['b',]), frozenset([])])
features = {(N-1,1): frozenset(['o',]),
            #(1,N-1): frozenset(['o',]),
            (3,2): frozenset(['h',]),
            (4,4): frozenset(['b',])}
heights = {(2,2): 1, (2,4): -1}
blur = 1
robot_nodes, real_robot_nodes = construct_nodes(l, N, label_set, features, heights, blur)

#-------------
#-------------

initial_node = (l, l, 'N')
initial_label = frozenset([])
home_xy = [(N, N)]
home_states = set([((2*n_x-1)*l, (2*n_x-1)*l, d) for (n_x, n_y) in home_xy for d in ['N', 'S', 'E', 'W']])

U = [tuple('FR'), tuple('BK'), tuple('TR'), tuple('TL')]
C = [3.0, 6.0, 5.0, 5.0]
P_FR = [1, 8, 1, 1]
P_BK = [2, 6, 2, 1]
P_TR = [1, 8, 1]
P_TL = [1, 8, 1]
P = [P_FR, P_BK, P_TR, P_TL]
robot_edges = construct_edges(robot_nodes, l, U, C, P)

visualize_world_paths(l, N, robot_nodes, [initial_node,initial_node,], [initial_label,initial_label], [], [0,0], 'test_map')

#-------------
print '---------- Construct robot mdp ----------'
robot_mdp = Motion_MDP(robot_nodes, robot_edges, U, initial_node, initial_label, home_states)
#robot_mdp.verify()
print '---------- Construct real robot mdp ----------'
real_mdp = Motion_MDP(real_robot_nodes, robot_edges, U, initial_node, initial_label, home_states)
t2 = time.time()
print '------------------------------'
print 'MDP done, time: %s' %str(t2-t0)


# ----
radius = 6
decay = 0.05
slope = [radians(-60), radians(45)]
robot_sensor = sensor(real_mdp, radius, decay, slope)


# ----
print '------------------------------'
base = 'G F b'
order = 'G i h X U ! h b'
safe = 'G ! o'
task1 = '& %s & %s %s' %(base, order, safe)
task = '& G F b G F h'
print 'Formula received: %s' %str(task)
dra = Dra(task)
t3 = time.time()

print 'DRA done, time: %s' %str(t3-t2)

#----
print '------------------------------'
gamma = [0.4, 0.9] # gamma_o, gamma_r
prod_dra = Product_Dra(mdp=robot_mdp, dra=dra, gamma=gamma)
#prod_dra.dotify()
t41 = time.time()
print 'Product DRA constructed, time: %s' %str(t41-t3)
#prod_dra.verify()
#----
prod_dra.init_dirichlet()
t42 = time.time()
print 'Compute init_dirichlet done, time: %s' %str(t42-t41)
#----
prod_dra.compute_init_mean_sigma()
t43 = time.time()
print 'Compute init_mean_sigma done, time: %s' %str(t43-t42)
#prod_dra.verify()
#----
prod_dra.compute_S_f()
t44 = time.time()
print 'Compute MEC done, time: %s' %str(t44-t43)

# #------
t5 = time.time()
total_T = 20
X, L, U, M, PX = execution_with_sensing(prod_dra, robot_sensor, total_T)
print '------------------------------'
print 'Planning and execution for %d steps, time: %s' %(total_T, str(t5-t44))
print 'Trajectory:', X
print 'Trace:', L
print 'Segment', M
print 'Action', U

# print '------------------------------'
# print 'Planning and execution for %d steps, time: %s' %(total_T, str(t5-t44))

