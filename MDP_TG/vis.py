import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

import numpy as np
import scipy.stats as stats
import random

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


def visualize_world_paths(l, Ns, robot_nodes, X, L, U, M, name=None):
    #----visualize simulated runs----
    k_colors = plt.get_cmap('copper')
    figure = plt.figure()
    ax = figure.add_subplot(1,1,1)
    #----- draw the workspace
    for node, prop in robot_nodes.iteritems():
        label = prop[0]
        if frozenset(['b']) in label.keys():
            text = '$b$'                
            color = 'yellow'
            density = label[frozenset(['b'])]
        elif frozenset(['o']) in label.keys():
            text = '$o$'
            color = 'red'
            density = label[frozenset(['o'])]
        elif frozenset(['h',]) in label.keys():
            text = '$h$'
            color = '#00bfff'
            density = label[frozenset(['h'])]
        elif frozenset(['w',]) in label.keys():
            text = '$w$'
            color = '#ADD8E6'
            density = label[frozenset(['w'])]
        else:
            text = None
            color = 'white'
        ht = prop[1]
        if color == 'white':
            if ((ht>=0.1*l) or (ht<=-0.1*l)):
                MIN_H, MAX_H = [-2*l, 1.5*l]
                c = (MAX_H - ht)/(MAX_H - MIN_H)
                color = k_colors(c)
        rec = matplotlib.patches.Rectangle((node[0]-l, node[1]-l),
                                           l*2, l*2,
                                           fill = True,
                                           facecolor = color,
                                           edgecolor = 'black',
                                           linewidth = 1,
                                           alpha =0.8)
        ax.add_patch(rec)
        if text:
            ax.text(node[0], node[1], r'%s' %text, fontsize = 10, fontweight = 'bold')
    K = len(X)
    #print 'K: %s' %K
    for k in xrange(0, K):
        l_wid = 2
        if M[k] == 0:
            Ecolor = 'blue'
        if M[k] == 1:
            Ecolor = 'magenta'
        if M[k] == 2:
            Ecolor = 'black'
        if M[k] > 2:
            Ecolor = 'green'
            l_wid = 10
        #----
        if (k<= K-2):
            line = matplotlib.lines.Line2D([X[k][0],X[k+1][0]],
                                           [X[k][1],X[k+1][1]],
                                           linestyle='-',
                                           linewidth=l_wid,
                                           color=Ecolor)
            ax.add_line(line)
            xl = X[k][0]
            yl = X[k][1]
            dl = X[k][2]
            if dl == 'N':
                car=[(xl-0.4*l,yl-0.4*l), (xl-0.4*l,yl+0.4*l), (xl, yl+0.8*l), (xl+0.4*l, yl+0.4*l), (xl+0.4*l,yl-0.4*l)]
            if dl == 'E':
                car=[(xl-0.4*l,yl+0.4*l), (xl+0.4*l,yl+0.4*l), (xl+0.8*l, yl), (xl+0.4*l, yl-0.4*l), (xl-0.4*l,yl-0.4*l)]
            if dl == 'S':
                car=[(xl+0.4*l,yl+0.4*l), (xl+0.4*l,yl-0.4*l), (xl, yl-0.8*l), (xl-0.4*l, yl-0.4*l), (xl-0.4*l,yl+0.4*l)]
            if dl == 'W':
                car=[(xl+0.4*l,yl-0.4*l), (xl-0.4*l,yl-0.4*l), (xl-0.8*l, yl), (xl-0.4*l, yl+0.4*l), (xl+0.4*l,yl+0.4*l)]
            polygon = Polygon(car, facecolor='black', edgecolor='black', lw=0.5, alpha=0.7, zorder = 4)
            ax.add_patch(polygon)
    ax.set_aspect('equal')
    ax.set_xlim(0, Ns*2*l)
    ax.set_ylim(0, Ns*2*l)
    ax.set_xlabel(r'$x(m)$')
    ax.set_ylabel(r'$y(m)$')    
    if name:
        plt.savefig('%s.pdf' %name, bbox_inches='tight')
        print "%s.pdf saved" %name


            
def analyze_events(MM, LL):
    #----analyze the results of all runs----
    #prefix failure/success, suffix failure/success 
    N = len(LL)
    T = len(LL[0])
    failure_count = 0.0
    prefix_failure_count = 0.0
    suffix_failure_count = 0.0
    prefix_suc_count = 0.0
    suffix_suc_count = 0.0
    for k, L in enumerate(LL):
        M = MM[k]
        Failure = False
        Prefix_suc = False
        Suffix_suc = False
        ip_count = 0
        for n,l in enumerate(L):
            m = M[n]
            if (m==2) and (not Suffix_suc):
                Failure = True
                failure_count += 1
                if not Prefix_suc:
                    prefix_failure_count += 1
                else:
                    suffix_failure_count += 1
                break
            if (m == 1) and (not Prefix_suc):
                Prefix_suc = True
                prefix_suc_count += 1
            if m == 10:
                ip_count += 1
            if ip_count >= 1:
                Suffix_suc = True
                suffix_suc_count += 1
                break
    print 'Analyze done'
    print 'Total %s simulations: %s failure (%s) (%s prefix, %s suffix), %s prefix successful (%s); %s suffix successful (%s)' %(str(N), str(failure_count), str(failure_count/N), str(prefix_failure_count), str(suffix_failure_count), str(prefix_suc_count), str(prefix_suc_count/N), str(suffix_suc_count), str(suffix_suc_count/prefix_suc_count))

    

def visualize_state_action_dynamic(motion_mdp, WS_d, WS_node_dict, x, l, u, m):
    #----visualize dynamic workspace and robot motion at each step----
    #----with action name and possible post states----
    figure = plt.figure()
    ax = figure.add_subplot(1,1,1)
    #----- draw the workspace
    xl = x[0]
    yl = x[1]
    dl = x[2]
    if m == 0:
        Ecolor = 'green'
    if m == 1:
        Ecolor = 'magenta'
    if m == 2:
        Ecolor = 'black'
    if m > 2:
        Ecolor = 'magenta'        
    if dl == 'N':
        car=[(xl-0.2,yl-0.2), (xl-0.2,yl+0.2), (xl, yl+0.4), (xl+0.2, yl+0.2), (xl+0.2,yl-0.2)]
    if dl == 'E':
        car=[(xl-0.2,yl+0.2), (xl+0.2,yl+0.2), (xl+0.4, yl), (xl+0.2, yl-0.2), (xl-0.2,yl-0.2)]
    if dl == 'S':
        car=[(xl+0.2,yl+0.2), (xl+0.2,yl-0.2), (xl, yl-0.4), (xl-0.2, yl-0.2), (xl-0.2,yl+0.2)]
    if dl == 'W':
        car=[(xl+0.2,yl-0.2), (xl-0.2,yl-0.2), (xl-0.4, yl), (xl-0.2, yl+0.2), (xl+0.2,yl+0.2)]                
    polygon = Polygon(car, fill = True, facecolor=Ecolor, edgecolor=Ecolor, lw=5, zorder=2)
    ax.add_patch(polygon)
    #
    actstr = r''
    for s in u:
        actstr += s
    ax.text(xl, yl+0.5, r'$%s$' %str(actstr), fontsize = 13, fontweight = 'bold', color='red')
    # plot shadow
    t_x_list = []
    for t_x in motion_mdp.successors(x):
        prop = motion_mdp[x][t_x]['prop']
        if u in prop.keys():
            t_x_list.append((t_x, prop[u][0]))
    #
    for new_x in t_x_list:
        xl = new_x[0][0]
        yl = new_x[0][1]
        dl = new_x[0][2]
        if dl == 'N':
            car=[(xl-0.2,yl-0.2), (xl-0.2,yl+0.2), (xl, yl+0.4), (xl+0.2, yl+0.2), (xl+0.2,yl-0.2)]
        elif dl == 'E':
            car=[(xl-0.2,yl+0.2), (xl+0.2,yl+0.2), (xl+0.4, yl), (xl+0.2, yl-0.2), (xl-0.2,yl-0.2)]
        elif dl == 'S':
            car=[(xl+0.2,yl+0.2), (xl+0.2,yl-0.2), (xl, yl-0.4), (xl-0.2, yl-0.2), (xl-0.2,yl+0.2)]
        elif dl == 'W':
            car=[(xl+0.2,yl-0.2), (xl-0.2,yl-0.2), (xl-0.4, yl), (xl-0.2, yl+0.2), (xl+0.2,yl+0.2)]
        polygon = Polygon(car, fill = True, facecolor='grey', edgecolor='grey', lw=5, zorder = 1)
        ax.add_patch(polygon)
        prob = new_x[1]
        ax.text(xl, yl, r'$%s$' %str(prob), fontsize = 10, fontweight = 'bold', color='red')                
    #----------------
    for node, prop in WS_node_dict.iteritems():
        if node != (x[0], x[1]):
            S = []
            P = []
            for s, p in prop.iteritems():
                S.append(s)
                P.append(p)
            rdn = random.random()
            pc = 0
            for k, p in enumerate(P):
                pc += p
                if pc> rdn:
                    break
            current_s = S[k]
        if node == (x[0], x[1]):
            current_s = l
        #------
        if current_s == frozenset(['base1','base']):
            text = '$base1$'
            color = 'yellow'
        elif current_s == frozenset(['base2','base']):
            text = '$base2$'
            color = 'yellow'
        elif current_s == frozenset(['base3','base']):
            text = '$base3$'
            color = 'yellow'
        elif current_s == frozenset(['obstacle','low']):
            text = '$Obs$'
            color = '#ff8000'
        elif current_s == frozenset(['obstacle','top']):
            text = '$Obs$'
            color = 'red'
        elif current_s == frozenset(['supply',]):
            text = '$Sply$'
            color = '#0000ff'
        else:
            text = None
            color = 'white'
        rec = matplotlib.patches.Rectangle((node[0]-WS_d, node[1]-WS_d),
                                               WS_d*2, WS_d*2,
                                               fill = True,
                                               facecolor = color,
                                               edgecolor = 'black',
                                               linewidth = 1,
                                               ls = '--',
                                               alpha =0.8)
        ax.add_patch(rec)
        if text:
            ax.text(node[0]-0.7, node[1], r'%s' %text, fontsize = 10, fontweight = 'bold')        
    ax.set_aspect('equal')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_xlabel(r'$x(m)$')
    ax.set_ylabel(r'$y(m)$')    
    return figure    

 
def run_movie(motion_mdp, WS_d, WS_node_dict, X, L, U, M):
    #----save simulation plot at each time step----
    #----for movie compilation, see mkmovie.sh for details----
    DPI = 500
    i = 0
    figure1 = visualize_state_action_dynamic(motion_mdp, WS_d, WS_node_dict, X[0], L[0], U[0], M[0])
    figure1.savefig('movie/frame%s.png' %i, dpi=DPI)
    plt.close()
    i += 1
    for k in xrange(1,len(X)):
        figure2 = visualize_state_dynamic(motion_mdp, WS_d, WS_node_dict, X[k], L[k], U[k], M[k])
        figure2.savefig('movie/frame%s.png' %i, dpi=DPI)
        plt.close()        
        i+=1
        figure1 = visualize_state_action_dynamic(motion_mdp, WS_d, WS_node_dict, X[k], L[k], U[k], M[k])
        figure1.savefig('movie/frame%s.png' %i, dpi=DPI)
        plt.close()                
        i +=1        
        plt.close()


def compute_suffix_mean_cost(UU, MM, COST):
    # record mean total cost of accepting cyclic path
    mean_cost_U = []
    for k, U in enumerate(UU):
        M = MM[k]
        c = 0.0
        t = 0.0
        for j,u in enumerate(U):
            if M[j] == 1:
                c += COST[u]
                t += 1.0
            elif M[j] == 10:
                mean_cost_U.append(c/t)
                c = 0.0
                t = 0.0
    print 'total number of cyclic path:%d' %len(mean_cost_U)
    mean_cost = sum(mean_cost_U)/len(mean_cost_U)
    return mean_cost
