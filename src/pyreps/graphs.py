import networkx as nx


def matching_graph(n):
    """Makes a matching graph since a complete graph.

    Args:
        n (int): A integer that to do the complete graph.

    Returns:
        networkx.classes.graph.Graph: The matching graph.

    Raises:
        NetworkXError: If n is a negative number.

    """
    k_n = nx.complete_graph(n)
    G = nx.Graph()
    for i in k_n.edges():
        G.add_node(i)
    w = []
    for i in k_n.edges():
        for j in k_n.edges():
            if ((j[0] not in i) and (j[1] not in i) and ((i, j) not in w) and ((j, i) not in w)):
                w.append((i, j))
                G.add_edge(i, j)
    return G
