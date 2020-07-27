from collections import defaultdict
import graph_handle

STR_PADDING_SIZE = 30
COUNT_STR_PADDING_SIZE = 10
SCORES_FILE_NAME = 'scores_genomes.txt'
ADDITIONAL_FILE_NAME = 'additional_data.txt'


def parse_args(args):
    num_of_args = len(args) - 1

    if num_of_args < 3:
        raise ValueError('One or more arguments are missing')

    scoring_matrix = read_scoring_matrix(args[1])

    sequences = {}
    for i in range(2, num_of_args + 1):
        sequence_id, sequence = read_seq_file(args[i])
        sequences[sequence_id] = sequence

    return sequences, scoring_matrix


def build_sequences_dict(sequences, K):
    mapped_sequences = {}
    for sequence_id, sequence in sequences.items():
        sequence_dict = map_sequence(sequence, K)
        mapped_sequences[sequence_id] = sequence_dict

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
            if score >= T:
                neighbor = list(neighbor)
                neighbor[pos] = char
                neighbor = "".join(neighbor)
                find_neighbors_rec(kmer, neighbor, pos + 1, score, alphabet, neighbors, scoring_matrix, T)


def create_file_for_final_scores(msps_dict, scores_list):
    f = open(SCORES_FILE_NAME, "w")
    counter = 0
    for sequences in msps_dict.keys():
        seq1 = sequences[0]
        seq2 = sequences[1]
        f.write(seq1 + "\t" + seq2 + "\t" + str(scores_list[counter]) + "\n")
        counter += 1

    f.close()


def write_additional_data(seq1, seq2, hsps_count, msps_count):
    fp = open(ADDITIONAL_FILE_NAME, "a")
    fp.write(('Sequence1 id: {} Sequence2 id: {} HSPs found: {} MSPs found: {}\n'
              .format(seq1.ljust(STR_PADDING_SIZE), seq2.ljust(STR_PADDING_SIZE),
                      str(hsps_count).ljust(COUNT_STR_PADDING_SIZE), str(msps_count).ljust(COUNT_STR_PADDING_SIZE))))


def calculate_scores(msps):
    scores_list = []
    for pair in msps.values():
        g = graph_handle.create_graph(pair)
        path = graph_handle.find_path(g)
        path_score = graph_handle.compute_pairwise_score(path, g)
        scores_list.append(path_score)

    return scores_list

