from collections import defaultdict

# TODO: write a function that parses the file 'sub_mat.txt'
def parse_subtitution_matrix(scoring_matrix_path):
    pass


def read_seq_file(seq_file):
    seq_id = ''
    seq = ''
    with open(seq_file) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
            else:
                seq += line.strip()

        return seq_id, seq


def map_sequence(sequence, k):
    sequence_dict = defaultdict(list)
    for i in range(0, len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        sequence_dict[kmer].append(i)
    return sequence_dict


def align(seq1, seq2, scoring_matrix):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be with the same length")

    score = 0

    for i in range(len(seq1)):
        score += scoring_matrix[seq1[i], seq2[i]]

    return score


def find_neighbors(kmer, scoring_matrix, alphabet, T):
    neighbors = []
    max_score = align(kmer, kmer, scoring_matrix)

    if max_score >= T:
        find_neighbors_rec(kmer, kmer, 0, max_score, alphabet, neighbors, scoring_matrix, T)

    return neighbors


def find_neighbors_rec(kmer, neighbor, pos, curr_score, alphabet, neighbors, scoring_matrix, T):
    if len(kmer) == pos:
        neighbors.append(neighbor)
    else:
        for char in alphabet:
            score = curr_score - scoring_matrix[kmer[pos], kmer[pos]] + scoring_matrix[kmer[pos], char]
            if (score) >= T:
                neighbor = list(neighbor)
                neighbor[pos] = char
                neighbor = "".join(neighbor)
                find_neighbors_rec(kmer, neighbor, pos + 1, score, alphabet, neighbors, scoring_matrix, T)


