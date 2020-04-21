import math
import itertools
import networkx as nx

from pyreps.utilities import tuple_sorted


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


def clique_graph(g, cmax=math.inf):
    """Makes a clique graph since a matching graph.

    Args:
        n (int): A integer that to do the matching graph.

    Returns:
        networkx.classes.graph.Graph: The clique graph.

    Raises:
        NetworkXError: If n is a negative number.

    """
    ite = nx.find_cliques(g)
    cliques = []
    K = nx.Graph()
    while True:
        try:
            cli = next(ite)
            cliques.append(frozenset(cli))
            if len(cliques) > cmax:
                return None
        except StopIteration:
            break
    K.add_nodes_from(cliques)
    clique_pairs = itertools.combinations(cliques, 2)
    K.add_edges_from((c1, c2) for (c1, c2) in clique_pairs if c1 & c2)
    G1 = nx.Graph()
    for i in K.nodes():
        G1.add_node(tuple(sorted(i)))
    for i in K.edges():
        if (tuple(sorted(i[0])) < tuple(sorted(i[1]))):
            G1.add_edge(tuple(sorted(i[0])), tuple(sorted(i[1])))
        else:
            e = tuple_sorted((tuple(sorted(i[1])), tuple(sorted(i[0]))))
            G1.add_edge(*e)
    return G1
