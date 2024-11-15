import pickle
import networkx as nx
from matplotlib import pyplot as plt

def get_equilibrium_nodes(nodes, reactions):
    equilibrium_nodes = {}
    for node in nodes:
        if isinstance(node, int) and node in reactions:
            equilibrium_nodes[node] = str(reactions[node]['e'])
    return equilibrium_nodes

def custom_layout(graph, node_colors):
    pos = {}
    left_x = 0.1
    right_x = 0.9
    y_step = 1 / (len(graph.nodes()) + 1)
    y = 1 - y_step
    for node, color in zip(graph.nodes(), node_colors):
        if color == 'skyblue':
            pos[node] = (left_x, y)
        else:
            pos[node] = (right_x, y)
        y -= y_step
    return pos

r = pickle.load(open("../app/data/SCAN_reactions.p", "rb"))

g = pickle.load(open("SCAN_graph3.p", "rb"))

reactions = r[300][10]
graph = g[300][10]

nodes = [node_id for node_id in graph.nodes()]
equilibrium_nodes = get_equilibrium_nodes(nodes, reactions)

# Create a color map for nodes
node_colors = ['green' if isinstance(node, int) else 'skyblue' for node in nodes]

# Create labels for nodes
labels = {node: equilibrium_nodes[node] if node in equilibrium_nodes else node for node in nodes}

plt.figure(figsize=(12, 8))
pos = custom_layout(graph, node_colors)
nx.draw(
    graph,
    pos,
    with_labels=True,
    labels=labels,
    node_size=300,
    node_color=node_colors,
    font_size=8,
    edge_color="gray",
)
plt.title("Simplified graph, first 100 nodes")
plt.show()



r1 = OrderedDict([('CH3CH2OH', 1), ('CO', 3), ('H2', 2), ('O2', 3)])