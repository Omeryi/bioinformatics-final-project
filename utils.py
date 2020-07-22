from collections import defaultdict

def parse_args(args):
    num_of_args = len(args) - 1

    if (num_of_args < 3):
        raise ValueError('One or more arguments are missing')

    scoring_matrix = read_scoring_matrix(args[1])

    sequences = {}
    for i in range(2, num_of_args + 1):
        sequence_name, sequence = read_seq_file(args[i])
        sequences[sequence_name] = sequence

    return sequences, scoring_matrix


def build_sequences_dict(sequences, K):
    mapped_sequences = {}
    for sequence_name, sequence in sequences.items():
        sequence_dict = map_sequence(sequence, K)
        mapped_sequences[sequence_name] = sequence_dict

    return mapped_sequences

def read_scoring_matrix(path):
    scoring_matrix = {}

    with open(path) as f:
        chars = f.readline().strip().split()

        for line in f:
            ch1, *scores = line.strip().split()

            for i, score in enumerate(scores):
                scoring_matrix[(ch1, chars[i])] = int(score)

    return scoring_matrix


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
        try:
            score += scoring_matrix[seq1[i],seq2[i]]
        except IndexError:
            print("i", i)
            print("seq1:" , seq1)
            print("seq2:", seq2)
            print("seq1[i]:", seq1[i])
            print("seq2[i]:", seq2[i])
            print("##################")

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


