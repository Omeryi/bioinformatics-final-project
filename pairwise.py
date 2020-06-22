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
    k = 5
    T = 24
    db = 'VTSAQELRGCTVINGSLIINIRGGNNLAAELEANLG'
    query = 'AANIEGGNNAAA'
    matrix = substitution_matrices.load("BLOSUM62")
    print("hi")

    hsps = hsp.get_hsps(query, db, k, matrix, T)
    print(hsps[0])

    X = 20
    msp = hsp.extend_hsp(query, db, hsp, matrix, X)
    print(msp)