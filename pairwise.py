import sys
import utils
import hsp_handle as hsp
from Bio.Align import substitution_matrices


# this function return a list of the sequences, each sequence is represented by tuple (id, seq)
def get_sequences():
    sequences = []
    [sequences.append(utils.read_seq_file(seq)) for seq in sys.argv[1:]]
    return sequences


if __name__ == '__main__':

    sequences, scoring_matrix = utils.parse_args(sys.argv)
    mapped_sequences = utils.build_sequences_dict(sequences)

    runtime_data = {}
    alignments = {}

