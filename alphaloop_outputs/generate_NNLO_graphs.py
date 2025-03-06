import yaml
import os
from pprint import pformat


def sort_by_id(graph_name_and_graph):
    graph_name, graph = graph_name_and_graph
    strip_GL = graph_name.split('GL')[1]
    return int(strip_GL)


graphs = []

category = 'NNLO_graphs'


files = [os.path.join(category, file)
         for file in os.listdir(category) if file.endswith('.yaml')]

for file in files:
    with open(file, 'r') as f:
        graphs.append([os.path.basename(file)[:-5],
                      yaml.load(f, Loader=yaml.FullLoader)])

graphs.sort(key=sort_by_id)

# load categories
contents = None
with open(os.path.join(category, 'categories.py'), 'r') as f:
    contents = f.read()
exec(contents)

print("nf_graphs=", nf_graphs)
print("s_channel_singlet=", s_channel_singlet)
print("t_channel_singlet=", t_channel_singlet)
print("outer_nf_graphs=", outer_nf_graphs)
print("no_valid_lmb=", no_valid_lmb)

# Fixes to graph
# for (g_name, g) in graphs:
#    g['edges'] = {k+1: v for k,v in g['edges'].items()}
#    g['nodes'] = {k+1: v for k,v in g['nodes'].items()}
#    for e in g['edges'].values():
#        e['vertices'] = [v+1 for v in e['vertices']]
#    for i_n, n in g['nodes'].items():
#        n['edge_ids'] = [e+1 for e in n['edge_ids']]
#        n['indices'] = [
#                g['edges'][e_id]['indices'][0] if (len(g['edges'][e_id]['indices'])==1 or g['edges'][e_id]['vertices'][0] == i_n)
#                else g['edges'][e_id]['indices'][1]
#            for e_id in n['edge_ids']
#        ]

inputs = [
    (
        'NNLO_from_gL__outer_nf', [
            (g_name, g) for (g_name, g) in graphs if (
                g_name in outer_nf_graphs and
                # this is already tested in rust, but we can do it again
                g_name not in s_channel_singlet and
                g_name not in t_channel_singlet)
        ],
    ),
    (
        'NNLO_from_gL__nf__not_singlet',
        [
            (g_name, g) for (g_name, g) in graphs if (
                g_name in nf_graphs and
                g_name not in outer_nf_graphs and
                g_name not in s_channel_singlet and
                g_name not in t_channel_singlet)
        ],
    ), (
        'NNLO_from_gL__non_nf__not_singlet',
        [
            (g_name, g) for (g_name, g) in graphs if (
                g_name not in nf_graphs and
                g_name not in s_channel_singlet and
                g_name not in t_channel_singlet)
        ],
    ), (
        'NNLO_from_gL__singlet__valid_lmb',
        [
            (g_name, g) for (g_name, g) in graphs if (
                (g_name in s_channel_singlet or
                 g_name in t_channel_singlet) and
                g_name not in no_valid_lmb
            )
        ],
    ), (
        'NNLO_from_gL__singlet__no_valid_lmb',
        [
            (g_name, g) for (g_name, g) in graphs if (
                (g_name in s_channel_singlet or
                 g_name in t_channel_singlet) and
                g_name in no_valid_lmb
            )
        ],
    )
]

# Sanity check
for (category, selected_graphs) in inputs:
    if category in ['NNLO_from_gL__nf__not_singlet', 'NNLO_from_gL__non_nf__not_singlet']:
        if any(g_name in no_valid_lmb for (g_name, _g) in selected_graphs):
            print("Non-singlet diagrams should always have a valid LMB!")
            sys.exit(1)

    summed_selected_sg_name_list = sorted([sg_name for (
        cat, selected_graphs) in inputs for (sg_name, _sg) in selected_graphs])
    sg_name_list = sorted([gn for (gn, g) in graphs])
    if summed_selected_sg_name_list != sg_name_list:
        print("Category selections do not add up to the original list of graphs!")
        sys.exit(1)

for (category, selected_graphs) in inputs:
    with open(f'{category}.py', 'w') as f:
        f.write("graphs={}\ngraph_names=[{}]".format(
            pformat([g for g_name, g in selected_graphs]),
            ','.join('"%s"' % g_name for g_name, g in selected_graphs)
        ))
