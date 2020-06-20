from copy import copy

class HSP():
    """Represents an alignment with score >=T between a kmer in the query and a kmer in the db"""

    def __init__(self, query_start, query_end, db_start, db_end, score):
        self.query_start = query_start
        self.query_end = query_end
        self.db_start = db_start
        self.db_end = db_end
        self.score = score

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            is_query_eq = self.query_start == other.query_start and self.query_end == other.query_end
            is_db_eq = self.db_start == other.db_start and self.db_end == other.db_end
            return is_query_eq and is_db_eq and self.score == other.score
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str((self.query_start, self.query_end, self.db_start, self.db_end, self.score)))

    def __str__(self):
        return f'''
                Score: {self.score}
                Query: [{self.query_start}, {self.query_end}]
                DB: [{self.db_start}, {self.db_end}]
                '''

    def __repr__(self):
        return self.__str__()


def get_hsps(query, db_dict, k, scoring_matrix, T):
    hsps = []
    query_dict = build_db(query, k)

    for key, val in query_dict.items():
        neighbors = find_neighbors(key, scoring_matrix, ALPHABET, T)
        for neighbor in neighbors:
            if neighbor in db_dict.keys():
                score = align(key, neighbor, scoring_matrix)
                for query_start in val:
                    query_end = query_start + k - 1
                    for db_start in db_dict[neighbor]:
                        db_end = db_start + k - 1
                        hsp = HSP(query_start, query_end, db_start, db_end, score)
                        hsps.append(hsp)

    return hsps


def extend_left(query, db, hsp, scoring_matrix, X):
    """returns the left extension with the maximal score"""

    msp = copy(hsp)

    while (hsp.score > msp.score - X and hsp.query_start > 0 and hsp.db_start > 0):
        hsp.query_start = hsp.query_start - 1
        hsp.db_start = hsp.db_start - 1
        hsp.score = hsp.score + scoring_matrix[query[hsp.query_start], db[hsp.db_start]]
        if (hsp.score > msp.score):
            msp = copy(hsp)

    return msp


def extend_right(query, db, hsp, scoring_matrix, X):
    """returns the right extension with the maximal score"""

    msp = copy(hsp)

    while (hsp.score > msp.score - X and hsp.query_end < len(query) - 1 and hsp.db_end < len(db) - 1):
        hsp.query_end = hsp.query_end + 1
        hsp.db_end = hsp.db_end + 1
        hsp.score = hsp.score + scoring_matrix[query[hsp.query_end], db[hsp.db_end]]
        if (hsp.score > msp.score):
            msp = copy(hsp)

    return msp


def extend_hsp(query, db, hsp, scoring_matrix, X):
    msp_left = extend_left(query, db, hsp, scoring_matrix, X)
    msp = extend_right(query, db, msp_left, scoring_matrix, X)

    return msp


def find_msps(query, db, k, scoring_matrix, T, X):
    msps = set()

    db_dict = build_db(db, k)
    hsps = get_hsps(query, db_dict, k, scoring_matrix, T)

    for hsp in hsps:
        msp = extend_hsp(query, db, hsp, scoring_matrix, X)
        msps.add(msp)

    return msps

