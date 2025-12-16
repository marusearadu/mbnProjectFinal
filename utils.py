import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import networkx as nx

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
    F1 = 2*TP/(2*TP + FP + FN)
    print(f"True Positives: \t{TP} \nFalse Negatives: \t{FN} \nFalse Positives: \t{FP} \nPrecision: \t\t{PR:.2f} \nRecall: \t\t{RE:.2f} \nF1: \t\t{F1:.2f}")
    
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

def plot_predicted_vs_actual(model, data, factors, title=''):
    # Plot predicted vs actual for each variable
    fig, ax = plt.subplots(3, 2, figsize=(15, 10))
    ax_flat = ax.ravel()
    total_mse = 0
    for idx, var in enumerate(factors.columns):
        variables, mu, cov = model.predict(factors.drop(columns=[var]))
        x = data['time']
        # Calculate MSE
        mse = np.mean((factors[var] - mu.ravel())**2)
        total_mse += mse
        print(f'MSE for {var}: {mse}')
        max_val = factors[var].max()
        max_val = max_val + 0.7 * abs(max_val)
        min_val = factors[var].min()
        min_val = min_val - 0.7 * abs(max_val)
        ax_flat[idx].plot(x, factors[var], label='Actual')
        ax_flat[idx].plot(x, mu.ravel(), label='Predicted')
        # Plot 95% confidence interval
        std = np.sqrt(cov).ravel()
        ax_flat[idx].fill_between(x, (mu.ravel() - 1.96 * std), (mu.ravel() + 1.96 * std), color='gray', alpha=0.5, label='95% CI')
        ax_flat[idx].set_ylim([min_val, max_val])
        ax_flat[idx].set_title(f'Predicted vs Actual for {var} {title}')
        ax_flat[idx].set_xlabel('Time')
        ax_flat[idx].set_ylabel(var)
        ax_flat[idx].grid()
        ax_flat[idx].legend()
    plt.tight_layout()
    plt.show()
    print(f'Total MSE: {total_mse}')

def deloop(net):
    flag = True
    while flag:
        try:
            edges = nx.find_cycle(net)
            # Remove a random edge from the cycle
            edge_to_remove = random.choice(edges)
            net.remove_edge(*edge_to_remove)
        except nx.exception.NetworkXNoCycle:
            flag = False
            pass
    return net