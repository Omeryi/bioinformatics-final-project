import sys
import utils
import hsp_handle as HSP
import time
import graph_handle


if __name__ == '__main__':
    start_time = time.time()

    K = 12
    T = 58
    X = 30

    sequences, scoring_matrix = utils.parse_args(sys.argv)
    mapped_sequences = utils.build_sequences_dict(sequences, K)
    msps = HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X)

    score_list = []
    for pair in msps.values():
         g = graph_handle.creating_graph(pair)
         path = graph_handle.find_path(g)
         path_score = graph_handle.compute_pairwise_score(path, g)
         score_list.append(path_score)

    utils.creating_file_for_final_scores(msps, score_list)

    total_runtime = (time.time() - start_time)
    fp = open("Additional_Data.txt", "a")
    fp.write(('\nTotal runtime of the program: {} Seconds\n'.format(total_runtime)))
    print("--- %s seconds ---" % total_runtime)









