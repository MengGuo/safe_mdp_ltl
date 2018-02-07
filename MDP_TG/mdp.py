# -*- coding: utf-8 -*-

from math import sqrt
from networkx.classes.digraph import DiGraph
from networkx import strongly_connected_component_subgraphs

from dirichlet import est_mean_sigma

#----------------------------------

class Motion_MDP(DiGraph):
    #----construct probabilistic-labeled MDP----
    def __init__(self, node_dict, edge_dict, U, initial_node, initial_label, home_states):
        DiGraph.__init__(self, name='motion_mdp', init_state=initial_node, init_label=initial_label, home=home_states)
        for (n, prob_label) in node_dict.iteritems():
            self.add_node(n, label = prob_label[0], height = prob_label[1], act = set())
        print "-------Motion MDP Initialized-------"
        self.add_edges(edge_dict, U)
        print "%s states and %s edges" %(str(len(self.nodes())), str(len(self.edges())))
        self.add_mean_sigma()
        
    def add_edges(self, edge_dict, U):
        self.graph['U'] = set()
        for u in U:
            self.graph['U'].add(tuple(u))
        for edge, attri in edge_dict.iteritems():
            f_node = edge[0]
            u = edge[1]
            t_node = edge[2]
            if (f_node, t_node) in self.edges():
                prop = self[f_node][t_node]['prop']
                prop[tuple(u)] = [attri[0], attri[1], 0, 0] #k,c,p,sigma
            else:                
                prob_cost = dict()
                #--- prob, cost
                prob_cost[tuple(u)] = [attri[0], attri[1], 0, 0]
                self.add_edge(f_node, t_node, prop = prob_cost)
        #----
        for f_node in self.nodes():
            Us = set()
            for t_node in self.successors(f_node):
                prop = self[f_node][t_node]['prop']
                Us.update(set(prop.keys()))
            if Us:
                self.node[f_node]['act'] = Us.copy()
            else:
                print 'Isolated state'
        print "-------Motion MDP Constructed-------"

    def set_gamma(self, gamma):
        # gamma_o, gamma_r
        self.graph['gamma'] = gamma    

    def add_mean_sigma(self):
        for f_node in self.nodes():
            for u in self.node[f_node]['act']:
                alpha = []
                b = []
                for t_node in self.successors(f_node):
                    prop = self[f_node][t_node]['prop']
                    if u in prop.keys():
                        alpha.append(prop[u][0])
                        b.append(t_node)
                # print 'alpha, b', [alpha, b]
                mean_b, sigma_b = est_mean_sigma(alpha, b)
                # print 'mean_b, sigma_b', [mean_b, sigma_b]
                for t_node in self.successors(f_node):
                    prop = self[f_node][t_node]['prop']
                    if u in prop.keys():
                        prop[u][2] = mean_b[t_node]
                        prop[u][3] = sigma_b[t_node]
        print 'Add mean and sigma Done'

    def verify(self):
        print '----------to verify MDP outgoing prob----------'
        for f_node in self.nodes():
            for u in self.node[f_node]['act']:
                print '----------'
                print 'from node %s under action %s' %(str(f_node), str(u))
                sum_p = 0
                for t_node in self.successors(f_node):
                    prop = self[f_node][t_node]['prop']
                    if u in prop.keys():
                        print 'to node %s: prop %s' %(str(t_node), str(prop[u]))
                        sum_p += prop[u][2]
                print 'total sum_p', sum_p
                    
#--------------------------------
def find_MECs(mdp, Sneg):
    #----implementation of Alg.47 P866 of Baier08----
    print 'Remaining states size', len(Sneg)
    U = mdp.graph['U']
    A = dict()
    for s in Sneg:
        A[s] = mdp.node[s]['act'].copy()
        if not A[s]:
            print "Isolated state"
    MEC = set()
    MECnew = set()
    MECnew.add(frozenset(Sneg))
    #----
    k = 0
    while MEC != MECnew:
        print "<============iteration %s============>" %k
        k +=1
        MEC = MECnew
        MECnew = set()
        print "MEC size: %s" %len(MEC)
        print "MECnew size: %s" %len(MECnew)
        for T in MEC:
            R = set()
            T_temp = set(T)
            simple_digraph = DiGraph()
            for s_f in T_temp:
                if s_f not in simple_digraph:
                    simple_digraph.add_node(s_f)
                for s_t in mdp.successors(s_f):
                    if s_t in T_temp:
                        simple_digraph.add_edge(s_f,s_t)
            print "SubGraph of one MEC: %s states and %s edges" %(str(len(simple_digraph.nodes())), str(len(simple_digraph.edges())))
            Sccs = strongly_connected_component_subgraphs(simple_digraph)
            i = 0
            for Scc in Sccs:
                i += 1
                if (len(Scc.edges())>=1):
                    for s in Scc.nodes():
                        U_to_remove = set() 
                        for u in A[s]:
                            for t in mdp.successors(s):
                                if ((u  in mdp[s][t]['prop'].keys()) and (t not in Scc.nodes())):
                                    U_to_remove.add(u)
                        A[s].difference_update(U_to_remove)
                        if not A[s]:                            
                            R.add(s)
            while R:
                s = R.pop()
                T_temp.remove(s)
                for f in mdp.predecessors(s):
                    if f in T_temp:
                        A[f].difference_update(set(mdp[f][s]['prop'].keys()))
                        if not A[f]:
                            R.add(f)
            New_Sccs = strongly_connected_component_subgraphs(simple_digraph)
            j = 0
            for Scc in New_Sccs:
                j += 1
                if (len(Scc.edges()) >= 1):
                    common = set(Scc.nodes()).intersection(T_temp)
                    if common:
                        MECnew.add(frozenset(common))
    #---------------
    print 'Final MEC and MECnew size:', len(MEC)
    return MEC, A



def find_SCCs(mdp, Sneg):
    #----simply find strongly connected components----
    print 'Remaining states size', len(Sneg)
    SCC  = set()
    simple_digraph = DiGraph()
    A = dict()
    for s in mdp.nodes():
        A[s] = mdp.node[s]['act'].copy()
    for s_f in Sneg:
        if s_f not in simple_digraph:
            simple_digraph.add_node(s_f)
        for s_t in mdp.successors(s_f):
            if s_t in Sneg:
                simple_digraph.add_edge(s_f,s_t)
    print "SubGraph of one Sf: %s states and %s edges" %(str(len(simple_digraph.nodes())), str(len(simple_digraph.edges())))
    sccs = strongly_connected_component_subgraphs(simple_digraph)
    for scc in sccs:
        SCC.add(frozenset(scc.nodes()))    
    return SCC, A

            
        
        
            
        
    
    
                    

