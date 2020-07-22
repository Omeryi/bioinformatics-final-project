import sys
import utils
import hsp_handle as HSP
from Bio.Align import substitution_matrices


# this function return a list of the sequences, each sequence is represented by tuple (id, seq)
def get_sequences():
    sequences = []
    [sequences.append(utils.read_seq_file(seq)) for seq in sys.argv[1:]]
    return sequences


if __name__ == '__main__':

    K = 11
    T = 20
    X = 100


    sequences, scoring_matrix = utils.parse_args(sys.argv)
    mapped_sequences = utils.build_sequences_dict(sequences, K)


    check = (HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X))
    print(check)
    #msps = HSP.create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X)

    #print(msps)
    runtime_data = {}
    alignments = {}
    print("done")



