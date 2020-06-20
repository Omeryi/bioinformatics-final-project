from collections import defaultdict

# TODO: write a function that parses the file 'sub_mat.txt'
def parse_subtitution_matrix(scoring_matrix_path):
    pass


def build_db(db, k):
    db_dict = defaultdict(list)
    for i in range(0, len(db) - k + 1):
        kmer = db[i:i + k]
        db_dict[kmer].append(i)
    return db_dict