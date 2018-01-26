

'''
Utilities to build quickly initial map and real map. 
'''


def blur_feature(n_xy, blur, features):
    # to find blurred feature 
    n_x = n_xy[0]
    n_y = n_xy[1]
    for key, value in features:
        f_x, f_y = key[:]        
        if ((f_x-blur<= n_x <= f_x+blur) and (f_y-blur<= n_y <= f_y+blur)):
            return value
    return None


def blur_height(n_xy, blur, heights):
    # to find blurred height
    n_x = n_xy[0]
    n_y = n_xy[1]
    for key, value in heights:
        f_x, f_y = key[:]        
        if ((f_x-blur<= n_x <= f_x+blur) and (f_y-blur<= n_y <= f_y+blur)):
            return value
    return None    



def construct_node(label_set, features, heights, blur):
    node_dict = dict()
    real_node_dict = dict()
    for n_x in range(N):
        for n_y in range(N):
            n_xy = (n_x+1, n_y+1)
            node_xy = (l*(2*n_x+1), l*(2*n_y+1))
            # ---------- initial map
            feat = blur_feature(n_xy, blur, features)
            label = dict()
            for l in label_set:
                if l != frozenset([]):
                    if l in feat:
                        label[l] = 3
                    else:
                        label[l] = 1
                else:
                    if not feat:
                        label[l] = 3
                    else:
                        label[l] = 1
            node_dict[node_xy] = [label, ]
            # ---------- real map
            real_label = dict()
            for l in label_set:
                if (n_xy in features):
                    if (l == features[node_xy]):
                        real_label[l] = 10
                    else:
                        real_label[l] = 1
                else:
                    if (l != frozenset([])):
                        real_label[l] = 1
                    else:
                        real_label[l] = 10
            real_node_dict[node_xy] = [real_label, ]
            # ---------- initial height
            ht = blur_height(n_xy, blur, heights)
            if ht:
                node_dict[node_xy].append(ht)
            else:
                node_dict[node_xy].append(0)
            # ---------- real height
            if n_xy in heights:
                real_node_dict[node_xy].append(heights[n_xy])
            else:
                real_node_dict[node_xy].append(0)
    return node_dict, real_node_dict
