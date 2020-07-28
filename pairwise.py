import sys
import utils
import hsp_handle as HSP
import time


if __name__ == '__main__':

    start_time = time.time()

    K = 12
    T = 58
    X = 30

    sequences, scoring_matrix = utils.parse_args(sys.argv)

    # For efficiency purposes, all sequences are mapped in advance
    mapped_sequences = utils.build_sequences_dict(sequences, K)
    msps = HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X)
    scores_list = utils.calculate_scores(msps)

    utils.create_file_for_final_scores(msps, scores_list)

    total_runtime = (time.time() - start_time)

    fp = open(utils.ADDITIONAL_FILE_NAME, "a")
    fp.write(('\nTotal runtime of the program: {} Seconds\n'.format(total_runtime)))










