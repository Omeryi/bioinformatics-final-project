
import utils
import itertools
from copy import copy
ALPHABET = 'ATGCSWRYKMBVHDNU'
# NOF:
# I changes a few names here so it makes more sense (db and query was less relevant)
# the functions here may be changed so that the build of the dict will only happen once (it takes too long)


class HSP():
    """Represents an alignment with score >=T between a kmer in the seq1 and a kmer in the seq2"""

    def __init__(self, seq1_start, seq1_end, seq2_start, seq2_end, score):
        self.seq1_start = seq1_start
        self.seq1_end = seq1_end
        self.seq2_start = seq2_start
        self.seq2_end = seq2_end
        self.score = score

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            is_seq1_eq = self.seq1_start == other.seq1_start and self.seq1_end == other.seq1_end
            is_seq2_eq = self.seq2_start == other.seq2_start and self.seq2_end == other.seq2_end
            return is_seq1_eq and is_seq2_eq and self.score == other.score
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str((self.seq1_start, self.seq1_end, self.seq2_start, self.seq2_end, self.score)))

    def __str__(self):
        return f''' 
                Score: {self.score} 
                Query: [{self.seq1_start}, {self.seq1_end}] 
                DB: [{self.seq2_start}, {self.seq2_end}] 
                '''

    def __repr__(self):
        return self.__str__()



def get_hsps(seq1, seq2, K, scoring_matrix, T):
    hsps = []
    seq1_dict = utils.map_sequence(seq1, K)
    seq2_dict = utils.map_sequence(seq2, K)

    for key, val in seq1_dict.items():
        neighbors = utils.find_neighbors(key, scoring_matrix, ALPHABET, T)
        for neighbor in neighbors:
            if neighbor in seq2_dict.keys():
                score = utils.align(key, neighbor, scoring_matrix)
                for seq1_start in val:
                    seq1_end = seq1_start + K - 1
                    for seq2_start in seq2_dict[neighbor]:
                        seq2_end = seq2_start + K - 1
                        hsp = HSP(seq1_start, seq1_end, seq2_start, seq2_end, score)
                        hsps.append(hsp)

    return hsps



def extend_left(seq1, seq2, hsp, scoring_matrix, X):
    """returns the left extension with the maximal score"""

    msp = copy(hsp)

    while (hsp.score > msp.score - X and hsp.seq1_start > 0 and hsp.seq2_start > 0):
        hsp.seq1_start = hsp.seq1_start - 1
        hsp.seq2_start = hsp.seq2_start - 1
        hsp.score = hsp.score + scoring_matrix[seq1[hsp.seq1_start], seq2[hsp.seq2_start]]
        if (hsp.score > msp.score):
            msp = copy(hsp)

    return msp


def extend_right(seq1, seq2, hsp, scoring_matrix, X):
    """returns the right extension with the maximal score"""

    msp = copy(hsp)

    while (hsp.score > msp.score - X and hsp.seq1_end < len(seq1) - 1 and hsp.seq2_end < len(seq2) - 1):
        hsp.seq1_end = hsp.seq1_end + 1
        hsp.seq2_end = hsp.seq2_end + 1
        hsp.score = hsp.score + scoring_matrix[seq1[hsp.seq1_end], seq2[hsp.seq2_end]]
        if (hsp.score > msp.score):
            msp = copy(hsp)

    return msp


def extend_hsp(seq1, seq2, hsp, scoring_matrix, X):
    msp_left = extend_left(seq1, seq2, hsp, scoring_matrix, X)
    msp = extend_right(seq1, seq2, msp_left, scoring_matrix, X)

    return msp


def find_msps(seq1, seq2, k, scoring_matrix, T, X):
    msps = set()
    hsps = get_hsps(seq1, seq2, k, scoring_matrix, T)

# there is a need to add counter (according to the assignment)
    for hsp in hsps:
        msp = extend_hsp(seq1, seq2, hsp, scoring_matrix, X)
        msps.add(msp)

    return msps

def create_msps_dict(scoring_matrix, sequences, mapped_sequences, K, T, X):
    msps_dict = {}
    counter = 0

    for pair in itertools.combinations(sequences.items(), r = 2):

        seq1_id = pair[0][0]
        seq1_dict = pair[0][1]
        seq2_id = pair[1][0]
        seq2_dict = pair[1][1]

        msps_dict[(seq1_id, seq2_id)] = find_msps(seq1_dict, seq2_dict, K, scoring_matrix, T, X)

    return msps_dict





