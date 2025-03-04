import yaml
import os
from pprint import pformat
graphs = []

category = 'LO_graphs'


files = [os.path.join(category, file)
         for file in os.listdir(category) if file.endswith('.yaml')]

for file in files:
    with open(file, 'r') as f:
        graphs.append([os.path.basename(file)[:-5],yaml.load(f, Loader=yaml.FullLoader)])

# Fixes to graph
for (g_name, g) in graphs:
    g['edges'] = {k+1: v for k,v in g['edges'].items()}
    g['nodes'] = {k+1: v for k,v in g['nodes'].items()}
    for e in g['edges'].values():
        e['vertices'] = [v+1 for v in e['vertices']]
    for i_n, n in g['nodes'].items():
        n['edge_ids'] = [e+1 for e in n['edge_ids']]
        n['indices'] = [
                g['edges'][e_id]['indices'][0] if (len(g['edges'][e_id]['indices'])==1 or g['edges'][e_id]['vertices'][0] == i_n)
                else g['edges'][e_id]['indices'][1]
            for e_id in n['edge_ids']
        ]

with open(f'LO_from_gL.py', 'w') as f:
    f.write("graphs={}\ngraph_names=[{}]".format(
        pformat([g for g_name, g in graphs]),
        ','.join('"%s"'%g_name for g_name, g in graphs)
    ))
