import math
import warnings
import numpy as np
from scipy.optimize import minimize, LinearConstraint
import baltic as bt

# taken from treesort script reassortment_utils and translated from dendropy to baltic

# host_reassortments() written by me 

def sibling_distance(parent_node):
    return parent_node.children[0].length + parent_node.children[1].length

def top_1(tree, ignore_top_edges):
    
    import math
    
    edge_cutoff = math.inf
    edge_lengths = sorted([node.length for node in tree.Objects if node.length])
    top_percentile = int(round(len(edge_lengths) * (1.0 - ignore_top_edges / 100)))
    edge_cutoff = edge_lengths[top_percentile - 1]
    
    return(edge_cutoff)

def compute_rea_rate_simple(annotated_tree, evol_rate, ignore_top_edges=1, subtrees=False):

    if ignore_top_edges > 0:
        edge_cutoff = top_1(annotated_tree, ignore_top_edges)
    
    rea_events = 0
    for node in annotated_tree.Objects:
        if node.parent is None:
            continue  # Skip the root
        if node.length and node.length >= edge_cutoff:
            continue  # Skip the longest edges
        if node.traits.get('is_reassorted'):
            rea_annotation = node.traits.get('rea')
            is_uncertain = all([g_str.startswith('_') for g_str in rea_annotation.split('-')])
            rea_events += 0.5 if is_uncertain else 1
    
    tree_length = sum(node.length for node in annotated_tree.Objects if node.parent and node.length < edge_cutoff)
    
    rea_rate = rea_events / tree_length * evol_rate if tree_length > 0 else 0.0
    
    if subtrees == False:
        return rea_rate
    else:
        return rea_rate, rea_events, tree_length

def likelihood_binary(x, rea_events, edge_lengths):
    if x < 1e-10:
        return np.inf
    func = 0
    for i in range(len(rea_events)):
        if edge_lengths[i] > 0:
            if rea_events[i] > 0:
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        func -= np.log(1 - np.exp(-1 * x * edge_lengths[i]))
                    except Warning:
                        func += np.inf
            else:
                func -= (-1 * x * edge_lengths[i])
    return func

def compute_rea_rate_binary_mle(annotated_tree, evol_rate, ref_seg_len=1700, subtrees=False):
    
    rea_events = []
    edge_lengths = []
    processed_uncertain = set()
    
    for node in annotated_tree.Objects:
        if node.parent is None:
            continue  # Skip the root
        is_uncertain = False
        if node.traits.get('is_reassorted'):
            rea_annotation = node.traits.get('rea')
            is_uncertain = all([g_str.startswith('_') for g_str in rea_annotation.split('-')])
            if not is_uncertain:
                rea_events.append(1)
        else:
            rea_events.append(0)
        
        if is_uncertain:

            siblings = [child for child in node.parent.children if child != node]

            sibling = siblings[0] # the tree is binary, so there should only be one sibling

            if sibling not in processed_uncertain:
                rea_events.append(1)
                edge_length = sibling_distance(node.parent)
                processed_uncertain.add(node)
            else:
                continue
                
        edge_length = node.length
        
        if edge_length > 1e-7:
            edge_lengths.append(edge_length / evol_rate)
        elif rea_events[-1] > 0:
            edge_lengths.append((1 / ref_seg_len) / evol_rate)
        else:
            edge_lengths.append(0)
    
    est = compute_rea_rate_simple(annotated_tree, evol_rate, ignore_top_edges=1)
    np_est = np.array([est])
    linear_constraint = LinearConstraint([[1]], [0])
    num_est = minimize(likelihood_binary, np_est, args=(rea_events, edge_lengths), tol=1e-9,
                       constraints=[linear_constraint])
    
    if subtrees == False:
        return num_est.x[0] if num_est.success else None
    else:
        return num_est.x[0] if num_est.success else None, rea_events, edge_lengths

def host_reassortments(trait_trees, trait_order, clock_rate, ref_seg_len, ignore_top_edges=1):
    
    host_frequencies_simple = {trait: 0 for trait in trait_order}
    host_frequencies_mle = {trait: 0 for trait in trait_order}
    subtree_rates_mle = {trait: [] for trait in trait_order}
    subtree_rates_simple = {trait: [] for trait in trait_order}

    for trait, subtrees in trait_trees.items():
        
        total_tree_length = 0
        total_rea_flt = 0
        total_rea_list = []
        total_edge_lengths = []
        
        for _, subtree in subtrees:
            
            tree_length = 0
            rea_events_flt = 0
            edge_lengths = []
            rea_events_list = []
            
            # skipping canine h3n8 since i care about canine h3n2
            if trait == "Canine" and _ == "Equine":
                continue
                
            est1, rea_events_flt, tree_length = compute_rea_rate_simple(subtree, clock_rate, ignore_top_edges, subtrees=True)
            est2, rea_events_list, edge_lengths = compute_rea_rate_binary_mle(subtree, clock_rate, ref_seg_len, subtrees=True)

            subtree_rates_simple[trait].append(est1)
            subtree_rates_mle[trait].append(est2)
            
            total_tree_length += tree_length
            total_rea_flt += rea_events_flt
            total_rea_list.extend(rea_events_list)
            total_edge_lengths.extend(edge_lengths)

        # multiplying by clockrate cancels out units to be reassortments per year
        rate = 0 if total_tree_length == 0 else (total_rea_flt / total_tree_length) * (clock_rate)
        host_frequencies_simple[trait] = rate
        
        
        np_est = np.array(host_frequencies_simple[trait])

        linear_constraint = LinearConstraint([[1]], [0])
    
        num_est = minimize(likelihood_binary, np_est, args=(total_rea_list, total_edge_lengths), tol=1e-9,
                           constraints=[linear_constraint])
        
        if num_est.success:
            host_frequencies_mle[trait] = (num_est.x[0])
        else:
            host_frequencies_mle[trait] = None

    # these host trees have no reassortments 
    # remove = ['Camel', 'Mink', 'Seal', 'Feline']
    # for key in remove:
    #     host_frequencies_simple.pop(key, None)
    #     host_frequencies_mle.pop(key, None)
        
    return(host_frequencies_simple, host_frequencies_mle, subtree_rates_simple, subtree_rates_mle)