
import random

k = 4

# Load data
with open('data/scop_lookup.tsv') as file:
    lines = file.readlines()
sids, cls = list(zip(*(line.split() for line in lines)))

# Folds
def get_fold(c):
    return '.'.join(c.split('.')[:2])
folds = list(set([get_fold(c) for c in cls]))

# Split by folds
random.seed(42)
random.shuffle(folds)

n = len(folds)
chunk_sizes = [int(n/k)] * k
for i in range(n - sum(chunk_sizes)):
    chunk_sizes[i] += 1

splits = []
for i, size in enumerate(chunk_sizes):
    a = sum(chunk_sizes[:i])
    splits.append(folds[a:a+chunk_sizes[i]])

# Split input data
for i, split in enumerate(splits):
    l = []
    for sid, cl in zip(sids, cls):
        if get_fold(cl) in split:
            l.append(sid + ' ' + cl)

    with open(f'tmp/fold_split{i}.txt', 'w') as file:
        file.write('\n'.join(l))

