import re
import numpy as np

def parse_cigar(cigar_string):
    ref, query = 0, 0
    matches = []

    for cnt, action in re.findall('([0-9]*)([IDMP])', cigar_string):
        cnt = int(cnt)

        if action == 'D':
            ref += cnt
        elif action == 'I':
            query += cnt
        elif action == 'M':  # use only Perfect matches
            ref += cnt
            query += cnt
        elif action == 'P':
            matches += [(ref + i, query + i) for i in range(cnt)]
            ref += cnt
            query += cnt
        else:
            raise ValueError(f'Action {action}')

    return np.array(matches)

def mutual_information(p_ab):
    p_a = p_ab.sum(axis=1)
    p_b = p_ab.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        log_scores = np.log2(p_ab / (p_a[:, np.newaxis] * p_b))
        return np.sum(p_ab * log_scores, where=np.isfinite(log_scores))