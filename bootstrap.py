
'''
Utilities to build quickly initialize map and real map. 
'''

import random

def blur_feature(n_xy, blur, features):
    # to find blurred feature
    feat = set()
    n_x = n_xy[0]
    n_y = n_xy[1]
    for key, value in features.iteritems():
        f_x, f_y = key[:]
        if (key == n_xy):
            feat.add(value)
        if ((f_x-blur<= n_x <= f_x+blur) and (f_y-blur<= n_y <= f_y+blur)):
            # if random.random() >= 0.3:
            #     feat.add(value)
            feat.add(value)
    return feat


def blur_height(n_xy, blur, heights):
    # to find blurred height
    n_x = n_xy[0]
    n_y = n_xy[1]
    for key, value in heights.iteritems():
        f_x, f_y = key[:]
        if (key == n_xy):
            return value
        if ((f_x-blur<= n_x <= f_x+blur) and (f_y-blur<= n_y <= f_y+blur)):
            # if random.random() >= 0.2:
            #     return value
            # else:
            #     return 0
            return value
    return 0    


def construct_nodes(L, N, label_set, features, heights, blur):
    ws_nodes = dict()
    real_ws_nodes = dict()
    for n_x in range(N):
        for n_y in range(N):
            n_xy = (n_x+1, n_y+1)
            node_xy = (L*(2*n_x+1), L*(2*n_y+1))
            # ---------- initial map
            feat = blur_feature(n_xy, blur, features)
            label = dict()
            for l in label_set:
                if l != frozenset([]):
                    if (feat and (l in feat)):
                        label[l] = 3
                else:
                    if not feat:
                        label[l] = 3
                    else:
                        label[l] = 1
            ws_nodes[node_xy] = [label, ]
            # ---------- real map
            real_label = dict()
            for l in label_set:
                if (n_xy in features):
                    if (l == features[n_xy]):
                        real_label[l] = 10
                else:
                    if (l == frozenset([])):
                        real_label[l] = 10
            real_ws_nodes[node_xy] = [real_label, ]
            # ---------- initial height
            ht = blur_height(n_xy, blur, heights)
            if ht:
                ws_nodes[node_xy].append(ht)
            else:
                ws_nodes[node_xy].append(0)
            # ---------- real height
            if n_xy in heights:
                real_ws_nodes[node_xy].append(heights[n_xy])
            else:
                real_ws_nodes[node_xy].append(0)
    robot_nodes = dict()
    real_robot_nodes = dict()
    for loc, prop in ws_nodes.iteritems():
        for d in ['N', 'S', 'E', 'W']:
            robot_nodes[(loc[0], loc[1], d)] = prop
    for loc, prop in real_ws_nodes.iteritems():
        for d in ['N', 'S', 'E', 'W']:
            real_robot_nodes[(loc[0], loc[1], d)] = prop                    
    return robot_nodes, real_robot_nodes


def construct_nodes_from_graph(mdp):
    robot_nodes = dict()
    for n in mdp.nodes():
        robot_nodes[n] = [mdp.node[n]['label'], mdp.node[n]['height']]
    return robot_nodes


def construct_edges(robot_nodes, l, U, C, P):
    P_FR, P_BK, P_TR, P_TL = P[:]
    robot_edges = dict()
    for fnode in robot_nodes.iterkeys():
        fx = fnode[0]
        fy = fnode[1]
        fd = fnode[2]
        # action FR
        u = U[0]
        c = C[0]
        if fd == 'N':
            t_nodes = [(fx-2*l, fy+2*l, fd), (fx, fy+2*l, fd), (fx+2*l, fy+2*l, fd), (fx, fy, fd)]
        if fd == 'S':
            t_nodes = [(fx-2*l, fy-2*l, fd), (fx, fy-2*l, fd), (fx+2*l, fy-2*l, fd), (fx, fy, fd)]
        if fd == 'E':
            t_nodes = [(fx+2*l, fy-2*l, fd), (fx+2*l, fy, fd), (fx+2*l, fy+2*l, fd), (fx, fy, fd)]
        if fd == 'W':
            t_nodes = [(fx-2*l, fy-2*l, fd), (fx-2*l, fy, fd), (fx-2*l, fy+2*l, fd), (fx, fy, fd)]
        for k, tnode in enumerate(t_nodes):
            if tnode in robot_nodes.keys():
                robot_edges[(fnode, u, tnode)] = (P_FR[k], c)
        # action BK
        u = U[1]
        c = C[1]
        if fd == 'N':
            t_nodes = [(fx-2*l, fy-2*l, fd), (fx, fy-2*l, fd), (fx+2*l, fy-2*l, fd), (fx, fy, fd)]
        if fd == 'S':
            t_nodes = [(fx-2*l, fy+2*l, fd), (fx, fy+2*l, fd), (fx+2*l, fy+2*l, fd), (fx, fy, fd)]
        if fd == 'E':
            t_nodes = [(fx-2*l, fy-2*l, fd), (fx-2*l, fy, fd), (fx-2*l, fy+2*l, fd), (fx, fy, fd)]
        if fd == 'W':
            t_nodes = [(fx+2*l, fy-2*l, fd), (fx+2*l, fy, fd), (fx+2*l, fy+2*l, fd), (fx, fy, fd)]
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
        #print 'robot_edges', robot_edges
    return robot_edges
