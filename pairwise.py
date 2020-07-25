import sys
import utils
import hsp_handle as HSP
import time
import graph_handle
# from Bio.Align import substitution_matrices


# this function return a list of the sequences, each sequence is represented by tuple (id, seq)
def get_sequences():
    sequences = []
    [sequences.append(utils.read_seq_file(seq)) for seq in sys.argv[1:]]
    return sequences


if __name__ == '__main__':
    start_time = time.time()

    # K = 4
    # T = 20
    # X = 100
    K = 12
    T = 60
    X = 20

    sequences, scoring_matrix = utils.parse_args(sys.argv)
    mapped_sequences = utils.build_sequences_dict(sequences, K)

    check = HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X)
    # print(check)
    # msps = HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X)

    # print(msps)
    runtime_data = {}
    alignments = {}
    print("done")
    print("--- %s seconds ---" % (time.time() - start_time))
    print("....................................................................................................")

# CODE FOR CHECKING THE FILES A. B. C
    # print()
    # print("--- CHECK: GRAPH SECTION ---")
    # g = graph_handle.creating_graph(check[('A', 'C')])
    # path = graph_handle.find_path(g)
    # path_score = graph_handle.compute_pairwise_score(path, g)
    # print("the path is: " + str(path) + "\n the score of the path is: " + str(path_score))
    # graph_handle.drawing_graph(g)


# CODE FOR CHECKING TASK 5
    print()
    print("--- CHECK: GRAPH SECTION ---")
    for pair in check.values():
        g = graph_handle.creating_graph(pair)
        path = graph_handle.find_path(g)
        path_score = graph_handle.compute_pairwise_score(path, g)
        print("the path is: " + str(path) + "\n the score of the path is: " + str(path_score))


# CODE FOR CREATING A FILE FOR THE FINAL SCORES
    # score_list = []
    # for pair in check.values():
    #     g = graph_handle.creating_graph(pair)
    #     path = graph_handle.find_path(g)
    #     path_score = graph_handle.compute_pairwise_score(path, g)
    #     score_list.append(path_score)
    #
    # utils.creating_file_for_final_scores(check, score_list)





