from MDP_TG.mdp import Motion_MDP
from MDP_TG.dra import Dra, Product_Dra
from MDP_TG.lp import syn_full_plan
from MDP_TG.vis import visualize_run
from MDP_TG.sense import sensor

from boostrap import construct_node


from math import radians

import time

t0 = time.time()

#-------- construct model -------
l = 1.5
N = 4
label_set = set([frozenset(['o',]), frozenset(['h',]), frozenset(['b',]), frozenset([])])
features = {(4,1): frozenset(['o',]),
            (3,2): frozenset(['h',]),
            (4,4): frozenset(['b',])}
heights = {(2,2): 10, (2,4): -10}
blur = 1
ws_nodes, real_ws_nodes = construct_node(label_set, features, height, blur)

#-------------
robot_nodes = dict()
real_robot_nodes = dict()
for loc, prop in ws_nodes.iteritems():
    for d in ['N', 'S', 'E', 'W']:
        robot_nodes[(loc[0], loc[1], d)] = prop
for loc, prop in real_ws_nodes.iteritems():
    for d in ['N', 'S', 'E', 'W']:
        real_robot_nodes[(loc[0], loc[1], d)] = prop        
#-------------

initial_node = (l, l, 'N')
initial_label = frozenset([])
home_xy = [(4, 1)]
home_states = set([((2*n_x-1)*l, (2*n_x-1)*l, d) for (n_x, n_y) in home_xy for d in ['N', 'S', 'E', 'W']])

U = [tuple('FR'), tuple('BK'), tuple('TR'), tuple('TL')]
C = [3.0, 6.0, 5.0, 5.0]
P_FR = [1, 8, 1]
P_BK = [2, 6, 2]
P_TR = [1, 8, 1]
P_TL = [1, 8, 1]
#-------------
robot_edges = dict()
for fnode in robot_nodes.iterkeys():
    fx = fnode[0]
    fy = fnode[1]
    fd = fnode[2]
    # action FR
    u = U[0]
    c = C[0]
    if fd == 'N':
        t_nodes = [(fx-2*l, fy+2*l, fd), (fx, fy+2*l, fd), (fx+2*l, fy+2*l, fd)]
    if fd == 'S':
        t_nodes = [(fx-2*l, fy-2*l, fd), (fx, fy-2*l, fd), (fx+2*l, fy-2*l, fd)]
    if fd == 'E':
        t_nodes = [(fx+2*l, fy-2*l, fd), (fx+2*l, fy, fd), (fx+2*l, fy+2*l, fd)]
    if fd == 'W':
        t_nodes = [(fx-2*l, fy-2*l, fd), (fx-2*l, fy, fd), (fx-2*l, fy+2*l, fd)]
    for k, tnode in enumerate(t_nodes):
        if tnode in robot_nodes.keys():
            robot_edges[(fnode, u, tnode)] = (P_FR[k], c)
    # action BK
    u = U[1]
    c = C[1]
    if fd == 'N':
        t_nodes = [(fx-2*l, fy-2*l, fd), (fx, fy-2*l, fd), (fx+2*l, fy-2*l, fd)]
    if fd == 'S':
        t_nodes = [(fx-2*l, fy+2*l, fd), (fx, fy+2*l, fd), (fx+2*l, fy+2*l, fd)]
    if fd == 'E':
        t_nodes = [(fx-2*l, fy-2*l, fd), (fx-2*l, fy, fd), (fx-2*l, fy+2*l, fd)]
    if fd == 'W':
        t_nodes = [(fx+2*l, fy-2*l, fd), (fx+2*l, fy, fd), (fx+2*l, fy+2*l, fd)]                
    for k, tnode in enumerate(t_nodes):
        if tnode in robot_nodes.keys():
            robot_edges[(fnode, u, tnode)] = (P_BK[k], c)
    # action TR
    u = U[2]
    c = C[2]
    if fd == 'N':
        t_nodes = [(fx, fy, 'N'), (fx, fy, 'E'), (fx, fy, 'S')]
    if fd == 'S':
        t_nodes = [(fx, fy, 'S'), (fx, fy, 'W'), (fx, fy, 'N')]
    if fd == 'E':
        t_nodes = [(fx, fy, 'E'), (fx, fy, 'S'), (fx, fy, 'W')]
    if fd == 'W':
        t_nodes = [(fx, fy, 'W'), (fx, fy, 'N'), (fx, fy, 'E')]
    for k, tnode in enumerate(t_nodes):
        if tnode in robot_nodes.keys():
            robot_edges[(fnode, u, tnode)] = (P_TR[k], c)
    # action TL
    u = U[3]
    c = C[3]
    if fd == 'S':
        t_nodes = [(fx, fy, 'S'), (fx, fy, 'E'), (fx, fy, 'N')]
    if fd == 'N':
        t_nodes = [(fx, fy, 'N'), (fx, fy, 'W'), (fx, fy, 'S')]
    if fd == 'W':
        t_nodes = [(fx, fy, 'W'), (fx, fy, 'S'), (fx, fy, 'E')]
    if fd == 'E':
        t_nodes = [(fx, fy, 'E'), (fx, fy, 'N'), (fx, fy, 'W')]
    for k, tnode in enumerate(t_nodes):
        if tnode in robot_nodes.keys():
            robot_edges[(fnode, u, tnode)] = (P_TL[k], c)
    # action ST
    u = U[4]
    c = C[4]
    if fd == 'S':
        t_nodes = [(fx, fy, 'W'), (fx, fy, 'S'), (fx, fy, 'E')]
    if fd == 'N':
        t_nodes = [(fx, fy, 'W'), (fx, fy, 'N'), (fx, fy, 'E')]
    if fd == 'W':
        t_nodes = [(fx, fy, 'S'), (fx, fy, 'W'), (fx, fy, 'N')]
    if fd == 'E':
        t_nodes = [(fx, fy, 'N'), (fx, fy, 'E'), (fx, fy, 'S')]   
    for k, tnode in enumerate(t_nodes):
        if tnode in robot_nodes.keys():
            robot_edges[(fnode, u, tnode)] = (P_ST[k], c)                    
#----
robot_mdp = Motion_MDP(robot_nodes, robot_edges, U, initial_node, initial_label, home_states)
real_mdp = Motion_MDP(real_robot_nodes, robot_edges, U, initial_node, initial_label, home_states)
t2 = time.time()
print '------------------------------'
print 'MDP done, time: %s' %str(t2-t0)


# ----
radius = 6
decay = 0.1
slope = [radians(-60), radians(45)]
robot_sensor = sensor(real_mdp, radius, decay, slope)


# ----
print '------------------------------'
base = 'G F b'
order = 'G i h X U ! h b'
safe = 'G ! o'
task = '& %s & %s %s' %(base, order, safe)
print 'Formula received: %s' %str(task)
dra = Dra(task)
t3 = time.time()

print 'DRA done, time: %s' %str(t3-t2)

#----
print '------------------------------'
gamma = [0.7, 0.8]
prod_dra = Product_Dra(mdp=robot_mdp, dra=dra, gamma=gamma)
#prod_dra.dotify()
t41 = time.time()
print 'Product DRA constructed, time: %s' %str(t41-t3)
#----
prod_dra.compute_S_f()
t42 = time.time()
print 'Compute MEC done, time: %s' %str(t42-t41)
#----
prod_dra.init_dirichlet()
t43 = time.time()
print 'Compute init_dirichlet done, time: %s' %str(t43-t42)
#----
prod_dra.compute_init_mean_sigma()
t44 = time.time()
print 'Compute init_mean_sigma done, time: %s' %str(t44-t43)

#------
total_T = 2
X, L, U, M, PX = prod_dra.execution_with_sensing(sensor, total_T)
t5 = time.time()
print '------------------------------'
print 'Planning and execution for %d steps, time: %s' %(total_T, str(t5-t42))

