#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import re
import sys
import json
import argparse
import baltic as bt


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--tree', type=str, required=True, help='path to treesort .tre output')
parser.add_argument('--outdir', type=str, required=True, help='path to output node data json')
args = parser.parse_args()


# In[ ]:


''' 
this function preps the treesort output to be readable by baltic

1. converts nexus format to nwk
2. replaces commas with "-" (where the reassorted segments are inferred)
3. remove the single quotation marks around TS_NODE_####
4. replaces ? with _ (where there is an undetermined reassortment event)
   
'''

def prep(qc_input):
        
    with open(qc_input, 'r') as file:
        nexus = file.read()
        
    start_idx = nexus.find('(')
    modified = nexus[start_idx:]
    
    end_idx = modified.find('END;')
    modified = modified[:end_idx]
        
    # removing commas between segments
    modified = re.sub(r'&rea="([^"]+)"', lambda match: f'&rea="{match.group(1).replace(",", "-")}"', modified)
    
    # removing quotation marks around node names
    modified = re.sub(r"'(TS_NODE_\d+)'", r'\1', modified)
    
    # replacing ? with _ so baltic can read it in
    modified = modified.replace('?', '_')
    
    with open("results/output.nwk", "w") as output_file:
        output_file.write(modified.strip())
        
    mytree = bt.loadNewick('results/output.nwk', absoluteTime= False)
    
    return(mytree)


# In[ ]:


# adapted from jordan ort's code, translated from phylo.bio to baltic

# keep diergence number?

def reassortment_counter(mytree, output):
    
    rea_dict = {}
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

    for k in mytree.Objects:
        if k.traits["is_reassorted"]:
            # Extract and clean reassorted segments
            raw_segments = k.traits["rea"]
            segment_names = [seg.split("(")[0] for seg in raw_segments.split("-")]
            reassorted_segments = f"{len(segment_names)} ({', '.join(segment_names)})"

            if k.is_node():
                rea_dict[k.traits.get("label")] = {
                    "Reassorted": "True",
                    "Reassorted Segments": reassorted_segments
                }
            elif k.is_leaf():
                rea_dict[k.name] = {
                    "Reassorted": "True",
                    "Reassorted Segments": reassorted_segments
                }
        else:
            key = (k.traits["label"] if k.is_node() else k.name)
            rea_dict[key] = {"Reassorted": "False"}

    branch_dict = {}
    for k in mytree.Objects:
        if k.traits["is_reassorted"]:
            # Use the cleaned reassorted segments from `rea_dict`
            if k.is_node():
                branch_dict[k.traits.get("label")] = {
                    "labels": {'Reassorted Segments': rea_dict[k.traits.get("label")]['Reassorted Segments']}
                }
            elif k.is_leaf():
                branch_dict[k.name] = {
                    "labels": {'Reassorted Segments': rea_dict[k.name]['Reassorted Segments']}
                }


    out_dict = {'nodes': rea_dict, 'branches': branch_dict}

    with open(output, 'w') as f:
        json.dump(out_dict, f)


# In[ ]:


mytree = prep(args.tree)


# In[ ]:


'''

treesort will sometimes infer "uncertain" reassortment events
since it is a bifurcated tree, each internal node has only 2 direct children
if both children have uncertain reassortment events, randomly assign one child to be reassorted
the other child is stripped of the reassortment event
we then call reassortment_counter to generate rea.json
the tree file stays the same but the rea.json file does not have any uncertanties
therefore, the summary.json and summary tree will also not have any uncertanties

'''


def parse_rea_string(rea_str):
    return [seg for seg in rea_str.strip().split("-") if seg]

def rebuild_rea_string(segments):
    return "-".join(segments) if segments else None

for k in mytree.Objects:
    if k.is_node():
        children = k.children
    
        # only look at nodes whose children are both reassorted 
        # since that is the first requirment for a possible uncertain rea event
        if not all(child.traits.get("is_reassorted") == 1 for child in children):
            continue
        
        seg_lists = []
        
        for child in children:
            raw_rea = child.traits.get("rea", "")
            seg_lists.append(parse_rea_string(raw_rea))

        # Identify uncertain reassortment segments (start with "_") in both children
        segs0_uncertain = set(seg for seg in seg_lists[0] if seg.startswith("_"))
        segs1_uncertain = set(seg for seg in seg_lists[1] if seg.startswith("_"))
        
        shared_uncertain = segs0_uncertain & segs1_uncertain

        # this randomly assigns each uncertain segment to a random child 
        for seg in shared_uncertain:
            # print(seg)
            stripped = seg.lstrip("_")
            chosen = random.choice(children)
            other = [c for c in children if c is not chosen][0]
            # print("chosen: " + chosen.name if chosen.is_leaf() else "chosen: " + chosen.traits["label"])
            # print("other: " + other.name if other.is_leaf() else "other: " + other.traits["label"])
            
            # update chosen: replace _SEG(x) with SEG(x)
            chosen_rea = parse_rea_string(chosen.traits.get("rea", ""))
            chosen_rea.remove(seg)
            chosen_rea.append(stripped)
            chosen.traits["rea"] = rebuild_rea_string(chosen_rea)

            # update other: remove the uncertain segment
            other_rea = parse_rea_string(other.traits.get("rea", ""))
            other_rea.remove(seg)
            other.traits["rea"] = rebuild_rea_string(other_rea)
            # print(chosen_rea)
            # print(other_rea)
            # print("\n")
        
        for child in children:
            if child.traits.get("is_reassorted") == 1:
                rea_str = child.traits.get("rea", "")
                
                if not rea_str:
                    child.traits["is_reassorted"] = 0
                    child.traits.pop("rea", None)

                
reassortment_counter(mytree, args.outdir)

