import yaml
import os

graphs = []

category = 'NNLO_graphs/nf_graphs_no_valid_lmb'

files = [os.path.join(category, file)
         for file in os.listdir(category) if file.endswith('.yaml')]

for file in files:
    with open(file, 'r') as f:
        graphs.append(yaml.load(f, Loader=yaml.FullLoader))

print(graphs)
