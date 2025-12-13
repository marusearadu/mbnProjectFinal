import networkx as nx
import matplotlib.pyplot as plt

POS = {
    "CBF1":  (-1, 1.5),
    "ASH1":  (-1, -1.5),
    "SWI5":  (1, 0),
    "GAL4":  (-0.66, 0.5),
    "GAL80": (-0.66, -0.5),
}

REF_NET   = nx.DiGraph()
REF_EDGES = [
    # ("GAL80", "GAL4"),  # this is not actually a transcriptional regulation link, just a protein-protein regulation link, so we omit it
    ("GAL4", "SWI5"), 
    ("SWI5", "ASH1"), 
    ("SWI5", "GAL80"), 
    ("SWI5", "CBF1"), 
    ("ASH1", "CBF1"), 
    ("CBF1", "GAL4")
]
REF_NET.add_edges_from(REF_EDGES)

def highlight_edges(graph, reference):
    edge_colors = {}
    graph_res = nx.DiGraph()
    graph_res.add_nodes_from(graph.nodes())
    graph_res.add_edges_from(graph.edges())
    graph_res.add_edges_from(reference.edges())
    
    for edge in graph_res.edges():
        if reference.has_edge(*edge) and graph.has_edge(*edge):
            edge_colors[edge] = 'green' # true positive
        elif reference.has_edge(*edge):
            edge_colors[edge] = 'red'   # false negative
        else:
            edge_colors[edge] = 'blue'  # false positive
    
    TP = len(list(edge for edge in edge_colors if edge_colors[edge] == 'green'))
    FN = len(list(edge for edge in edge_colors if edge_colors[edge] == 'red'))
    FP = len(list(edge for edge in edge_colors if edge_colors[edge] == 'blue'))
    PR = TP / max(1, TP + FP)
    RE = TP / max(1, TP + FN)
    print(f"True Positives: \t{TP} \nFalse Negatives: \t{FN} \nFalse Positives: \t{FP} \nPrecision: \t\t{PR:.2f} \nRecall: \t\t{RE:.2f}")
    
    # Draw nodes
    fig, ax = plt.subplots(figsize=(12, 9))
    nx.draw_networkx_nodes(graph_res, POS, node_size=700, node_color='lightblue')
    nx.draw_networkx_labels(graph_res, POS, font_size=10, font_weight='bold')
    
    # Draw edges with curvature
    for edge, color in edge_colors.items():
        nx.draw_networkx_edges(
            graph_res, POS, 
            edgelist=[edge],
            edge_color=color,
            connectionstyle='arc3,rad=0.02',  # 0.2 gives anticlockwise curve
            arrows=True,
            arrowsize=20
        )
    
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], color = 'green', lw = 2, label = 'Common edges'),
        Line2D([0], [0], color = 'red',   lw = 2, label = 'Only in reference'),
        Line2D([0], [0], color = 'blue',  lw = 2, label = 'Only in graph')
    ]
    plt.legend(handles=legend_elements, loc='best')

    plt.show()