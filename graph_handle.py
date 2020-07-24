import networkx as nx
import hsp_handle


def find_path(graph):
    return nx.dag_longest_path(graph, weight='weight', default_weight=1, topo_order=None)


def compute_weight(score1, score2):
    return score1 + score2


def check_if_adding_edge(G, msp1, msp2):
    weight = compute_weight(msp1.score, msp2.score)
    G.add_edge(msp1, msp2, weight=weight)


def creating_graph(msps):
    msps_list = list(msps)
    G = nx.DiGraph()
    for i in range(len(msps_list)-1):
        for j in range(i+1, len(msps_list)):
            # if msp j starts after msp i
            if msps_list[i].seq1_end < msps_list[j].seq1_start and msps_list[i].seq2_end < msps_list[j].seq2_start:
                msp1 = msps_list[i]
                msp2 = msps_list[j]

            # if msp i starts after msp j
            elif msps_list[j].seq1_end < msps_list[i].seq1_start and msps_list[j].seq2_end < msps_list[i].seq2_start:
                msp1 = msps_list[j]
                msp2 = msps_list[i]

            else:
                continue
            check_if_adding_edge(G, msp1, msp2)
    return G


def compute_pairwise_score(path):
    total_score = 0
    for node in path:
        total_score += node.score
    return total_score


if __name__ == '__main__':
    hsp1 = hsp_handle.HSP(1, 3, 5, 6, 10)
    hsp2 = hsp_handle.HSP(7, 8, 20, 26, 5)
    g = creating_graph(set([hsp1, hsp2]))
    path = find_path(g)
    path_score = compute_pairwise_score(path)
    print("the path is: " + str(path) + "\n the score of the path is: " + str(path_score))

