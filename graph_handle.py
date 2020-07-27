import networkx as nx
import matplotlib.pyplot as plt


def find_path(graph):
    return nx.dag_longest_path(graph, weight='weight', default_weight=1, topo_order=None)


# The higher the scores of the 2 MSPs and the closer they are, the heavier their edge is
def compute_weight(score1, score2, distance_seq1, distance_seq2):
    return (score1 + score2) / (distance_seq1 + distance_seq2)


def calculate_distance(msp1, msp2):
    # msp1 always will start before msp2
    distance1 = msp2.seq1_start - msp1.seq1_end
    distance2 = msp2.seq2_start - msp1.seq2_end
    return distance1, distance2


def add_edge(G, msp1, msp2):
    # msp1 always will start before msp2
    distance1, distance2 = calculate_distance(msp1, msp2)
    weight = compute_weight(msp1.score, msp2.score, distance1, distance2)
    G.add_edge(msp1, msp2, weight=weight)


def create_graph(msps):
    msps_list = list(msps)
    G = nx.DiGraph()
    G.add_nodes_from(msps)
    for i in range(len(msps_list)-1):
        # msp1 always will start before msp2
        for j in range(i+1, len(msps_list)):
            # If msp j starts after msp i
            if msps_list[i].seq1_end < msps_list[j].seq1_start and msps_list[i].seq2_end < msps_list[j].seq2_start:
                msp1 = msps_list[i]
                msp2 = msps_list[j]

            # If msp i starts after msp j
            elif msps_list[j].seq1_end < msps_list[i].seq1_start and msps_list[j].seq2_end < msps_list[i].seq2_start:
                msp1 = msps_list[j]
                msp2 = msps_list[i]

            else:
                continue
            add_edge(G, msp1, msp2)

    return G


def compute_pairwise_score(path, G):
    total_score = sum(node.score for node in path)

    # Handling with the case in which it is preferred to take one high score MSP instead of a few concatenated MSPs
    msp_max_score = max(node.score for node in G.nodes())

    return max(total_score, msp_max_score)


# this function can be used for debug purposes
def draw_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size=500)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edges(graph, pos, edge_color='r', arrows=True)
    plt.show()

