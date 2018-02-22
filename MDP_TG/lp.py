# -*- coding: utf-8 -*-

from networkx.classes.digraph import DiGraph
from networkx import single_source_shortest_path
from gurobipy import *

from collections import defaultdict
import random


def syn_real_time_plan(prod_mdp, gamma, init_node):
    try:
        prod_mdp.v_return = syn_return_plan(prod_mdp)
        best_all_plan, segment = syn_full_plan(prod_mdp, gamma, init_node)
        return best_all_plan, segment
    except KeyboardInterrupt:
        print '-----interrupted plan synthesis!-----'
        return None, None

def syn_full_plan(prod_mdp, gamma, init_node, alpha=1):
    #----Optimal plan synthesis, total cost over plan prefix and suffix----
    print "==========[Optimal full plan synthesis start]=========="
    Plan = []
    segment = 'prefix'
    print 'number of S_fi: %d' %len(prod_mdp.Sf)
    for l, S_fi in enumerate(prod_mdp.Sf):
        print "---for one S_fi---"
        for k, MEC in enumerate(S_fi):
            sf = MEC[0]
            if init_node not in sf:
                print '----- init node in PREFIX ------'
                plan_prefix, prefix_cost, prefix_risk, y_in_sf, Sr, Sd = syn_plan_prefix(prod_mdp, MEC, gamma, init_node)
                print "Best plan prefix obtained, cost: %s, risk %s" %(str(prefix_cost), str(prefix_risk))
                plan_suffix, suffix_cost, suffix_risk = None, None, None
            else:
                print '----- init node in SUFFIX ------'
                segment = 'suffix'
                plan_suffix, suffix_cost, suffix_risk = syn_plan_suffix_new(prod_mdp, MEC, gamma, init_node)
                print "Best plan suffix obtained, cost: %s, risk %s" %(str(suffix_cost), str(suffix_risk))
                plan_prefix, prefix_cost, prefix_risk, y_in_sf, Sr, Sd = None, None, None, None, None, None
            if plan_prefix  or plan_suffix:
                plan = [[plan_prefix, prefix_cost, prefix_risk, y_in_sf], [plan_suffix, suffix_cost, suffix_risk], [MEC[0], MEC[1], Sr, Sd]]
                return plan, segment
    else:
        print "No valid plan found"
        return None, None

            

def syn_plan_prefix(prod_mdp, MEC, gamma, init_node):
    #----Synthesize optimal plan prefix to reach accepting MEC or SCC----
    #----with bounded risk and minimal expected total cost----
    print "===========[plan prefix synthesis starts]==========="
    gamma_o = gamma[0]
    gamma_r = gamma[1]
    Sf = MEC[0]
    ip = MEC[1]
    delta = 1.0
    print 'init_node', init_node
    print 'init node not in Sf', (init_node not in Sf)
    if init_node not in Sf:
        print '----- init node in PREFIX ------'
        path_init = single_source_shortest_path(prod_mdp, init_node)
        print 'Reachable from init size: %s' %len(path_init.keys())
        if not set(path_init.keys()).intersection(Sf):
            print "Initial node can not reach Sf"
            return None, None, None, None, None, None
        Sn = set(path_init.keys()).difference(Sf)
        print 'Reachable nodes from init and not in Sf, size: %s' %str(len(Sn))
        #----find bad states that can not reach MEC
        simple_digraph = DiGraph()
        simple_digraph.add_edges_from(((v,u) for u,v in prod_mdp.edges()))
        path = single_source_shortest_path(simple_digraph, random.sample(ip,1)[0])
        reachable_set = set(path.keys())
        print 'States that can reach Sf, size: %s' %str(len(reachable_set))
        Sd = Sn.difference(reachable_set)
        Sr = Sn.intersection(reachable_set)
        # #--------------
        print 'Sn size: %s; Sd inside size: %s; Sr inside size: %s' %(len(Sn),len(Sd), len(Sr))
        # ---------solve lp------------
        print '-----'
        print 'Gurobi starts now'
        print '-----'
        try:
            Y = defaultdict(float)
            model = Model('plan_prefix')
            # create variables
            for s in Sr:
                for u in prod_mdp.node[s]['act'].copy():
                    Y[(s,u)] = model.addVar(vtype=GRB.CONTINUOUS,lb=0, name='y[(%s,%s)]' %(s, u))
            model.update()
            print 'Variables added'
            # set objective
            obj = 0
            for s in Sr:
                for u in prod_mdp.node[s]['act']:
                    for t in prod_mdp.successors(s):
                        prop = prod_mdp[s][t]['prop'].copy()
                        if u in prop.keys():
                            xi = prop[u][3]
                            pe = prop[u][0]
                            ce = prop[u][2]
                            obj += Y[(s,u)]*pe*(ce-xi)
            model.setObjective(obj, GRB.MINIMIZE)
            print 'Objective function set'
            # add constraints
            #------------------------------
            y_to_sd = 0.0
            y_to_Sf = 0.0
            for s in Sr:
                for t in prod_mdp.successors(s):
                    prop = prod_mdp[s][t]['prop'].copy()
                    for u in prop.iterkeys():
                        sigma = prop[u][1]
                        pe = prop[u][0]
                        if t in Sd:
                            y_to_sd += Y[(s,u)]*pe*(1+sigma)
                        elif t in Sf:
                            y_to_Sf += Y[(s,u)]*pe*(1+sigma)
            # model.addConstr(y_to_Sf+y_to_sd >= delta, 'sum_out')
            model.addConstr(y_to_Sf >= gamma_o, 'risk')
            print 'Risk constraint added'
            #------------------------------
            y_to_return = 0.0
            p_to_return = 0.0
            for t in prod_mdp.successors(init_node):
                prop = prod_mdp[init_node][t]['prop'].copy()
                for u in prop.iterkeys():
                    sigma = prop[u][1]
                    v = prod_mdp.v_return[t[0]]
                    pe = prop[u][0]
                    y_to_return += Y[(init_node,u)]*pe*(v+sigma)
                    p_to_return += Y[(init_node,u)]
            model.addConstr(y_to_return >= gamma_r*p_to_return, 'safety')
            #model.addConstr(y_to_return >= gamma_r, 'safety')
            print 'Safety constraint added'
            #--------------------
            for t in Sr:
                node_y_in = 0.0
                node_y_out = 0.0
                for u in prod_mdp.node[t]['act']:   
                    node_y_out += Y[(t,u)]
                for f in prod_mdp.predecessors(t):
                    if f in Sr:
                        prop = prod_mdp[f][t]['prop'].copy()
                        for uf in prop.iterkeys():
                            #print 'f, uf, t, prop[uf]', [f, uf, t, prop[uf]]
                            node_y_in += Y[(f,uf)]*prop[uf][0]
                if t == init_node:
                    model.addConstr(node_y_out == 1.0 + node_y_in, 'init_node_flow_balance for %s' %str(t))
                    # model.addConstr(node_y_out == 1.0, 'init_node_out_flow for %s' %str(t))
                    # node_y_in = 0.0
                    # for f in prod_mdp.predecessors(t):
                    #     if (f in Sr) and (f != t):
                    #         prop = prod_mdp[f][t]['prop'].copy()
                    #         for uf in prop.iterkeys():
                    #             #print 'f, uf, t, prop[uf]', [f, uf, t, prop[uf]]
                    #             node_y_in += Y[(f,uf)]*prop[uf][0]
                    # model.addConstr(node_y_in == 0.0, 'init_node_in_flow for %s' %str(t))
                else:
                    model.addConstr(node_y_out == node_y_in, 'middle_node_flow_balance_for_%s' %str(t))
            print 'Initial node flow balanced'
            print 'Middle node flow balanced'
            model.update()
            #model.write("debug.lp")
            #----------------------
            # solve
            print '--optimization starts--'
            model.optimize()
            # print '--variables value--'
            # for v in model.getVars():
            #     print v.varName, v.x
            # print 'obj:', model.objVal
            #------------------------------
            # model.computeIIS()
            # print "The following constraint(s) cannot be satisfied:"
            # for c in model.getConstrs():
            #     if c.IISConstr:
            #         print c.constrName
            # compute plan prefix given the LP solution
            plan_prefix = dict()
            for s in Sr:
                norm = 0
                U = []
                P = []
                U_total = prod_mdp.node[s]['act'].copy()          
                for u in U_total:
                    norm += Y[(s,u)].X
                for u in U_total:
                    U.append(u)
                    if norm > 0.01:
                        P.append(Y[(s,u)].X/norm)
                    else:
                        P.append(1.0/len(U_total))
                plan_prefix[s] = [U, P]
            print "----Prefix plan generated----"
            print "----solution for current_node %s----" %str(init_node)
            Y_init_node = dict()
            for u in prod_mdp.node[init_node]['act'].copy():
                Y_init_node[u] =  Y[(init_node,u)].X
            print 'Y_init_node', Y_init_node
            print "----Plan for current_node %s----" %str(init_node)
            print plan_prefix[init_node]
            print "----Plan for current_node----"
            print plan_prefix[init_node]
            cost = model.objval            
            print "----Prefix cost computed"
            # compute the risk given the plan prefix
            risk = 0.0
            y_to_sd = 0.0
            y_to_Sf = 0.0
            for s in Sr:
                for t in prod_mdp.successors(s):
                    if t in Sd:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys(): 
                            pe = prop[u][0]
                            y_to_sd += Y[(s,u)].X*pe
                    elif t in Sf:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys(): 
                            pe = prop[u][0]
                            y_to_Sf += Y[(s,u)].X*pe
                            if Y[(s,u)].X > 0.0:
                                print '(s,u), Y[(s,u)].X, pe', [(s,u), Y[(s,u)].X, pe]
            if (y_to_sd+y_to_Sf) >0:
                risk = y_to_sd/(y_to_sd+y_to_Sf)
            print 'y_to_sd: %s; y_to_Sf: %s' %(y_to_sd, y_to_Sf)
            print "----Prefix risk computed: %s" %str(risk)
            # compute the input flow to the suffix
            y_in_Sf = dict()
            for s in Sn:
                for t in prod_mdp.successors(s):
                    if t in Sf:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys():
                            pe = prop[u][0]
                            if t not in y_in_Sf:
                                y_in_Sf[t] = Y[(s,u)].X*pe
                            else:
                                y_in_Sf[t] += Y[(s,u)].X*pe
            # normalize the input flow
            y_total = 0.0
            for s, y in y_in_Sf.iteritems():
                y_total += y
            print 'actual y_total: %s' %str(y_total)
            if y_total >0:
                for s in y_in_Sf.iterkeys():
                    y_in_Sf[s] = y_in_Sf[s]/y_total
                print "----Y in Sf computed and normalized"
            #print y_in_Sf                
            # compute safety
            y_to_return = []
            p_to_return = []
            for t in prod_mdp.successors(init_node):
                prop = prod_mdp[init_node][t]['prop'].copy()
                for u in prop.iterkeys():
                    sigma = prop[u][1]
                    v = prod_mdp.v_return[t[0]]
                    pe = prop[u][0]
                    y_to_return.append(Y[(init_node,u)].X*pe*(v+sigma))
                    p_to_return.append(Y[(init_node,u)].X)
            print "----Prefix safety computed: *%s* (%s) of *%s* (%s) = %.2f for given gamma_r %s" %(sum(y_to_return), str(y_to_return), sum(p_to_return), str(p_to_return), (sum(y_to_return)/sum(p_to_return)), gamma_r)
            return plan_prefix, cost, risk, y_in_Sf, Sr, Sd
        except GurobiError:
            print "Gurobi Error reported"
            return None, None, None, None, None, None
    else:
        print '----- init node not in PREFIX ------'
        return None, None, None, None, None, None


def syn_plan_suffix_new(prod_mdp, MEC, gamma, init_node):
    #----Synthesize optimal plan suffix to stay within the accepting MEC----
    #----with minimal expected total cost of accepting cyclic paths----
    print "===========[plan suffix synthesis starts]==========="
    gamma_o = gamma[0]
    gamma_r = gamma[1]
    Sf = MEC[0]
    ip = MEC[1]
    act = MEC[2].copy()
    y_in_Sf = dict()
    for s in Sf:
        #y_in_Sf[s] = 1.0/len(Sf)
        if s == init_node:
            y_in_Sf[s] = 1.0
        else:
            y_in_Sf[s] = 0.0
    if init_node in Sf:
        print '---------- init_node in SUFFIX ----------'
        if init_node in ip:
            print 'init_node in Ip'
        else:
            print 'init_node NOT in Ip'
        paths = single_source_shortest_path(prod_mdp, init_node)
        Sn = set(paths.keys()).intersection(Sf)
        print 'Sf size: %s' %len(Sf)
        print 'Sn as reachable Sf size: %s' %len(Sn)
        print 'Ip size: %s' %len(ip)
        print 'Ip and Sf intersection size: %s'%len(Sn.intersection(ip))
        # ---------solve lp------------
        print '------'
        print 'Gurobi starts now'
        print '------'
        try:
            Y = defaultdict(float)
            model = Model('plan_suffix')
            # create variables
            for s in Sn:
                for u in act[s]:
                    Y[(s,u)] = model.addVar(vtype=GRB.CONTINUOUS,lb=0, name='y[(%s, %s)]' %(s, u))
            model.update()
            print 'Variables added'
            # set objective
            obj = 0
            for s in Sn:
                for u in act[s]:
                    for t in prod_mdp.successors(s):
                        prop = prod_mdp[s][t]['prop'].copy()
                        if (u in prop.keys()):
                            pe = prop[u][0]
                            xi = prop[u][3]
                            ce = prop[u][2]
                            obj += Y[(s,u)]*pe*(ce - xi)
            model.setObjective(obj, GRB.MINIMIZE)
            print 'Objective added'
            # add constraints
            #--------------------
            for s in Sn:
                constr3 = 0
                constr4 = 0
                for u in act[s]:
                    constr3 += Y[(s,u)]
                for f in prod_mdp.predecessors(s):
                    if (f in Sn) and (s not in ip):
                        prop = prod_mdp[f][s]['prop'].copy()
                        for uf in act[f]:
                            if uf in prop.keys():
                                #print 'f, uf, s, prop[uf]', [f, uf, s, prop[uf]]
                                constr4 += Y[(f,uf)]*prop[uf][0]
                            else:
                                constr4 += Y[(f,uf)]*0.00
                    if (f in Sn) and (s in ip) and (f != s):
                        prop = prod_mdp[f][s]['prop'].copy()
                        for uf in act[f]:
                            if uf in prop.keys():
                                #print 'f, uf, s, prop[uf]', [f, uf, s, prop[uf]]
                                constr4 += Y[(f,uf)]*prop[uf][0]
                            else:
                                constr4 += Y[(f,uf)]*0.00
                if (s in y_in_Sf.keys()) and (s not in ip):
                    #print 's, in Sf, not in ip, y_in_Sf[s]', [s, y_in_Sf[s]]
                    model.addConstr(constr3 == constr4 + y_in_Sf[s], 'balance_with_y_in')
                if (s in y_in_Sf.keys()) and (s in ip):
                    #print 's, in Sf, in ip', [s, y_in_Sf[s]]
                    model.addConstr(constr3 == y_in_Sf[s], 'balance_with_y_in')
                if (s not in y_in_Sf.keys()) and (s not in ip):
                    #print 's, not in Sf, not in ip', s
                    model.addConstr(constr3 == constr4, 'balance')
            print 'Balance condition added'
            print 'Initial Sf condition added'
            #------------------------------
            y_to_ip = 0.0
            y_out = 0.0
            for s in Sn:
                for t in prod_mdp.successors(s):
                    if t not in Sn:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys():
                            if u in act[s]:
                                pe = prop[u][0]
                                y_out += Y[(s,u)]*pe
                    elif t in ip:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys():
                            if u in act[s]:
                                pe = prop[u][0]
                                y_to_ip += Y[(s,u)]*pe        
            #model.addConstr(y_to_ip >= gamma_o*(y_to_ip+y_out), 'risk')
            model.addConstr(y_to_ip >= gamma_o, 'risk')
            print 'Task constraint added'
            #------------------------------
            y_to_return = 0.0
            p_to_return = 0.0
            for t in prod_mdp.successors(init_node):
                prop = prod_mdp[init_node][t]['prop'].copy()
                for u in prop.iterkeys():
                    if u in act[init_node]:
                        sigma = prop[u][1]
                        v = prod_mdp.v_return[t[0]]
                        pe = prop[u][0]
                        y_to_return += Y[(init_node,u)]*pe*(v+sigma)
                        p_to_return += Y[(init_node,u)]
            #model.addConstr(y_to_return >= gamma_r*p_to_return, 'safety')
            model.addConstr(y_to_return >= gamma_r, 'safety')
            print 'Safety constraint added'
            #------------------------------
            # solve
            model.optimize()
            # print '--variables value--'
            # for v in model.getVars():
            #     print v.varName, v.x
            # print 'obj:', model.objVal
            #------------------------------
            # compute optimal plan suffix given the LP solution
            plan_suffix = dict()
            for s in Sn:
                norm = 0
                U = []
                P = []
                for u in act[s]:
                    norm += Y[(s,u)].X
                for u in act[s]:
                    U.append(u)
                    if norm > 0.01:
                        P.append(Y[(s,u)].X/norm)
                    else:
                        P.append(1.0/len(act[s])) 
                plan_suffix[s] = [U, P]
            # for key, value in plan_suffix.iteritems():
            #     print 'key, value', [key, value]
            print "----Suffix plan generated----"
            if init_node in ip:
                print 'init_node in Ip'
            else:
                print 'init_node NOT in Ip'
            print "----solution for current_node %s----" %str(init_node)
            Y_init_node_out = dict()
            for u in act[init_node]:
                Y_init_node_out[u] =  Y[(init_node,u)].X
            print 'Y_init_node_out', Y_init_node_out
            Y_init_node_in = dict()
            for f in prod_mdp.predecessors(init_node):
                if f in Sf:
                    prop = prod_mdp[f][init_node]['prop'].copy()
                    for uf in prop:
                        if uf in act[f]:
                            if prop[uf][0]*Y[(f,uf)].X>0:
                                Y_init_node_in[(f,uf)] = (prop[uf][0], Y[(f,uf)].X, prop[uf][0]*Y[(f,uf)].X)
            print 'Y_init_node_in', Y_init_node_in
            print "----Plan for current_node %s----" %str(init_node)
            print plan_suffix[init_node]
            cost = model.objval
            print "----Suffix cost computed"
            # compute task constraint given the plan suffix
            print '---now computing task constraint----'
            y_to_ip = 0.0
            for s in Sn:
                for t in prod_mdp.successors(s):
                    if t in ip:
                        prop = prod_mdp[s][t]['prop'].copy()
                        for u in prop.iterkeys():
                            if u in act[s]:
                                # print '[s, u, t, Y[(s,u)].X]', [s, u, t, Y[(s,u)].X]
                                pe = prop[u][0]
                                y_to_ip += Y[(s,u)].X*pe
                                if Y[(s,u)].X > 0.0:
                                    print '(s,u,t), Y[(s,u)].X, pe', [(s,u,t), Y[(s,u)].X, pe]
            print 'Task constraint: %s' %y_to_ip
            risk = y_to_ip
            # compute safety
            y_to_return = 0.0
            p_to_return = 0.0
            for t in prod_mdp.successors(init_node):
                prop = prod_mdp[init_node][t]['prop'].copy()
                for u in prop.iterkeys():
                    if u in act[init_node]:
                        sigma = prop[u][1]
                        v = prod_mdp.v_return[t[0]]
                        pe = prop[u][0]
                        y_to_return += Y[(init_node,u)].X*pe*(v+sigma)
                        p_to_return += Y[(init_node,u)].X
                        # print '[init_node, u, t, y_to_return, p_to_return]', [init_node, u, t, Y[(init_node,u)].X*pe*(v+sigma), Y[(init_node,u)].X]
            print "----Suffix safety computed: %s of %s = %.2f for given gamma_r %s" %(str(y_to_return), str(p_to_return), (y_to_return/p_to_return), gamma_r)
            return plan_suffix, cost, risk
        except GurobiError:
            print "Gurobi Error reported"
            return None, None, None
    else:
        print 'No initial node specified in suffix'
        return None, None, None, None, None, None

        
def syn_plan_suffix(prod_mdp, MEC, gamma, init_node):
    #----Synthesize optimal plan suffix to stay within the accepting MEC----
    #----with minimal expected total cost of accepting cyclic paths----
    print "===========[plan suffix synthesis starts]==========="
    gamma_o = gamma[0]
    gamma_r = gamma[1]
    sf = MEC[0]
    ip = MEC[1]
    act = MEC[2].copy()
    delta = 1.0
    if init_node in sf:
        print '---------- init_node in SUFFIX ----------'
        paths = single_source_shortest_path(prod_mdp, init_node)
        Sn = set(paths.keys()).intersection(sf)
        print 'Sf size: %s' %len(sf)
        print 'Sn as reachable sf size: %s' %len(Sn)
        print 'Ip size: %s' %len(ip)
        print 'Ip and sf intersection size: %s'%len(Sn.intersection(ip))
        Sr = Sn.difference(ip)
        # ---------solve lp------------
        print '------'
        print 'Gurobi starts now'
        print '------'
        try:
            if init_node not in ip:
                print '----- init_node NOT inside ip -----'
                Y = defaultdict(float)
                model = Model('plan_suffix')
                # create variables
                for s in Sr:
                    for u in act[s]:
                        Y[(s,u)] = model.addVar(vtype=GRB.CONTINUOUS,lb=0, name='y[(%s, %s)]' %(s, u))
                model.update()
                print 'Variables added'
                # set objective
                obj = 0
                for s in Sr:
                    for u in act[s]:
                        for t in prod_mdp.successors(s):
                            prop = prod_mdp[s][t]['prop'].copy()
                            if (u in prop.keys()):
                                pe = prop[u][0]
                                xi = prop[u][3]
                                ce = prop[u][2]
                                obj += Y[(s,u)]*pe*(ce - xi)
                model.setObjective(obj, GRB.MINIMIZE)
                print 'Objective added'
                # add constraints
                #--------------------
                for s in Sr:
                    print '--------------------'
                    node_y_out = 0.0
                    node_y_in = 0.0
                    for u in act[s]:
                        node_y_out += Y[(s,u)]
                    for f in prod_mdp.predecessors(s):
                        if ((f in Sr) and (f !=s)):
                            prop = prod_mdp[f][s]['prop'].copy()
                            for uf in prop.keys():
                                print 'f, uf, s, prop[uf]', [f, uf, s, prop[uf]]
                                if uf in act[f]:
                                    node_y_in += Y[(f,uf)]*prop[uf][0]
                    if s == init_node:
                        #model.addConstr(node_y_out == node_y_in + 1.0, 'init_flow_balance')
                        model.addConstr(node_y_out == 1.0, 'init_flow_balance')
                    else:
                        model.addConstr(node_y_out == node_y_in, 'flow_balance')
                print 'Balance condition added'
                print 'Initial sf condition added'
                #------------------------------
                y_to_ip = 0.0
                for s in Sr:
                    for t in prod_mdp.successors(s):
                        if t in ip:
                            prop = prod_mdp[s][t]['prop'].copy()
                            for u in prop.iterkeys():
                                if u in act[s]:
                                    pe = prop[u][0]
                                    y_to_ip += Y[(s,u)]*pe
                #model.addConstr(y_to_ip >= 0.01, 'task constraint')
                print 'Task constraint added'
                #------------------------------
                y_to_return = 0.0
                p_to_return = 0.0
                for t in prod_mdp.successors(init_node):
                    prop = prod_mdp[init_node][t]['prop'].copy()
                    for u in prop.iterkeys():
                        if u in act[init_node]:
                            sigma = prop[u][1]
                            v = prod_mdp.v_return[t[0]]
                            pe = prop[u][0]
                            y_to_return += Y[(init_node,u)]*pe*(v+sigma)
                            p_to_return += Y[(init_node,u)]
                #model.addConstr(y_to_return >= gamma_r*p_to_return, 'safety')
                #model.addConstr(y_to_return >= 0.0, 'safety')
                #model.addConstr(p_to_return == 1.0, 'safety')
                print 'Safety constraint added'
                #------------------------------
                # solve
                model.optimize()
                # print '--variables value--'
                for v in model.getVars():
                    print v.varName, v.x
                # print 'obj:', model.objVal
                #------------------------------
                # compute optimal plan suffix given the LP solution
                plan_suffix = dict()
                for s in Sr:
                    norm = 0
                    U = []
                    P = []
                    for u in act[s]:
                        norm += Y[(s,u)].X
                    for u in act[s]:
                        U.append(u)
                        if norm > 0.01:
                            P.append(Y[(s,u)].X/norm)
                        else:
                            P.append(1.0/len(act[s])) 
                    plan_suffix[s] = [U, P]
                for key, value in plan_suffix.iteritems():
                    print 'key, value', [key, value]                    
                cost = model.objval
                print "----Suffix cost computed"
                # compute task constraint given the plan suffix
                print '---now computing task constraint----'
                y_to_ip = 0.0
                for s in Sr:
                    for t in prod_mdp.successors(s):
                        if t in ip:
                            prop = prod_mdp[s][t]['prop'].copy()
                            for u in prop.iterkeys():
                                if u in act[s]:
                                    print '[s, u, t, Y[(s,u)].X]', [s, u, t, Y[(s,u)].X]
                                    pe = prop[u][0]
                                    y_to_ip += Y[(s,u)].X*pe
                print 'Task constraint: %s' %y_to_ip
                risk = y_to_ip
                # compute safety
                y_to_return = 0.0
                p_to_return = 0.0
                for t in prod_mdp.successors(init_node):
                    prop = prod_mdp[init_node][t]['prop'].copy()
                    for u in prop.iterkeys():
                        if u in act[init_node]:
                            sigma = prop[u][1]
                            v = prod_mdp.v_return[t[0]]
                            pe = prop[u][0]
                            y_to_return += Y[(init_node,u)].X*pe*(v+sigma)
                            p_to_return += Y[(init_node,u)].X
                            print '[init_node, u, t, y_to_return, p_to_return]', [init_node, u, t, Y[(init_node,u)].X*pe*(v+sigma), Y[(init_node,u)].X]
                print "----Suffix safety computed: %s of %s for given gamma_r %s" %(str(y_to_return), str(p_to_return), gamma_r)
                return plan_suffix, cost, risk
        except GurobiError:
            print "Gurobi Error reported"
            return None, None, None
    else:
        print 'No initial node specified in suffix'
        return None, None, None, None, None, None

        
def act_by_plan(prod_mdp, best_plan, segment, prod_state):
    #----choose the randomized action by the optimal policy---- 
    #recall that {best_plan = [plan_prefix, prefix_cost, prefix_risk, y_in_sf],
    #[plan_suffix, suffix_cost, suffix_risk], [MEC[0], MEC[1], Sr, Sd], plan_bad]}
    plan_prefix = best_plan[0][0]
    plan_suffix = best_plan[1][0]
    if (segment == 'prefix'):
        #print 'In prefix'
        U = plan_prefix[prod_state][0]
        P = plan_prefix[prod_state][1]
        rdn = random.random()
        pc = 0
        for k, p in enumerate(P):
            pc += p
            if pc>rdn:
                break
        #print 'action chosen: %s' %str(U[k])
        return U[k], 0
    elif (segment == 'suffix'):
        #print 'In suffix'
        U = plan_suffix[prod_state][0]
        P = plan_suffix[prod_state][1]
        rdn = random.random()
        pc = 0
        for k, p in enumerate(P):
            pc += p
            if pc>rdn:
                break
        #print 'action chosen: %s' %str(U[k])
        if prod_state in best_plan[2][1]:
            return U[k], 10
        else:
            return U[k], 1
    

def rd_act_by_plan(prod_mdp, best_plan, prod_state, I):
    #----choose the randomized action by the optimal policy----
    #optimal prefix with round-robin as the plan suffix
    #recall that {best_plan = [plan_prefix, prefix_cost, prefix_risk, y_in_sf],
    #[plan_suffix, suffix_cost, suffix_risk], [MEC[0], MEC[1], Sr, Sd], plan_bad]}
    plan_prefix = best_plan[0][0]
    plan_suffix = best_plan[1][0]
    plan_bad = best_plan[3]
    if (prod_state in plan_prefix):
        #print 'In prefix'
        U = plan_prefix[prod_state][0]
        P = plan_prefix[prod_state][1]
        rdn = random.random()
        pc = 0
        for k, p in enumerate(P):
            pc += p
            if pc>rdn:
                break
        #print 'action chosen: %s' %str(U[k])
        return U[k], 0, I
    elif (prod_state in plan_suffix):
        #print 'In suffix'
        U = plan_suffix[prod_state][0]
        #print 'action chosen: %s' %str(U[k])
        if I[prod_state] <= (len(U)-1):
            u = U[I[prod_state]]
            I[prod_state] += 1
        if I[prod_state] == len(U):
            I[prod_state] = 0
            u = U[I[prod_state]]
        if prod_state in best_plan[2][1]:
            return u, 10, I
        else:
            return u, 1, I
    elif (prod_state in plan_bad):
        #print 'In bad states'        
        U = plan_bad[prod_state][0]
        P = plan_bad[prod_state][1]
        rdn = random.random()
        pc = 0
        for k, p in enumerate(P):
            pc += p
            if pc>rdn:
                break
        #print 'action chosen: %s' %str(U[k])
        return U[k], 2, I

        
def syn_return_plan(prod_mdp):
    #----Synthesize optimal return plan, given set of home state S_h----
    print "===========[return plan synthesis starts]==========="
    # ---------solve lp------------
    print '-----'
    print 'Gurobi starts now'
    print '-----'
    try:
        mdp = prod_mdp.graph['mdp']
        s_h = mdp.graph['home']
        V = defaultdict(float)
        model = Model('return_plan')
        # create variables
        for s in mdp.nodes():
            V[s] = model.addVar(vtype=GRB.CONTINUOUS,lb=0, name='v[%s]' %str(s))
        model.update()
        print 'Variables added'
        # --------------------
        # set objective
        obj = 0
        for s in mdp.nodes():
            obj += V[s]
        model.setObjective(obj, GRB.MINIMIZE)
        print 'Objective function set'
        # add constraints
        #------------------------------        
        #--------------------
        for s in mdp.nodes():
            #print 's now', s
            if s in s_h:
                model.addConstr(V[s] == 1.0, 'goal_node_value')
            else:
                for u in mdp.node[s]['act']:
                    sigma = 0
                    suc_V = 0
                    for t in mdp.successors(s):
                        prop = mdp[s][t]['prop']
                        if u in prop.keys():
                            # print 's,u,t,prop[u]',[s,u,t,prop[u]]
                            sigma += prop[u][2]*prop[u][3]
                            suc_V += prop[u][2]*V[t]
                    model.addConstr(V[s] >= sigma + suc_V, 'bellman_balance')
        print 'Goal node value added'
        print 'Bellman balanced'
        #----------------------
        #----------------------
        # solve
        print '--optimization starts--'
        model.optimize()
        # print '--variables value--'
        # for v in model.getVars():
        #     print v.varName, v.x
        print 'obj:', model.objVal
        #------------------------------
        print '--optimization done--'
        V_return = dict()
        for s in mdp.nodes():
            V_return[s] = V[s].X
        print '--value function returned--'
        # print 'V_return', V_return
        return V_return
    except GurobiError:
        print "Gurobi Error reported"
        return None
