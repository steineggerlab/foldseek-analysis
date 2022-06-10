#! /usr/bin/env python3
"""
    ./create_vqvae_training_data.py data/pdbs 270 0 2
"""

import os.path
import numpy as np
import random
import sys

import extract_pdb_features
import util

data_dir = os.path.join(os.path.dirname(__file__), 'data/')


feature_cache = {}  # path: (features, mask)
def encoder_features(pdb_path, virt_cb):
    """
    Calculate 3D descriptors for each residue of a PDB file.
    """
    feat = feature_cache.get(pdb_path, None)
    if feat is not None:
        return feat

    coords, valid_mask = extract_pdb_features.get_coords_from_pdb(pdb_path, full_backbone=True)
    coords = extract_pdb_features.move_CB(coords, virt_cb=virt_cb)
    partner_idx = extract_pdb_features.find_nearest_residues(coords, valid_mask)
    features, valid_mask2 = extract_pdb_features.calc_angles_forloop(coords, partner_idx, valid_mask)

    seq_dist = (partner_idx - np.arange(len(partner_idx)))[:, np.newaxis]
    log_dist = np.sign(seq_dist) * np.log(np.abs(seq_dist) + 1)

    vae_features = np.hstack([features, log_dist])
    feature_cache[pdb_path] = vae_features, valid_mask2

    return vae_features, valid_mask2


def align_features(pdb_dir, virtual_center, sid1, sid2, cigar_string):
    """
    Return aligned descriptors for a given alignment between two PDBs.
    """

    idx_1, idx_2 = util.parse_cigar(cigar_string).T

    feat1, mask1 = encoder_features(os.path.join(pdb_dir, sid1), virtual_center)
    feat2, mask2 = encoder_features(os.path.join(pdb_dir, sid2), virtual_center)

    valid_mask = mask1[idx_1] & mask2[idx_2]
    idx_1 = idx_1[valid_mask]
    idx_2 = idx_2[valid_mask]

    x = np.vstack([feat1[idx_1], feat2[idx_2]])
    y = np.vstack([feat2[idx_2], feat1[idx_1]])
    return x, y  # (n x 10, n x 10)


if __name__ == '__main__':
    pdb_dir = sys.argv[1]
    pairfile = sys.argv[2]
    a = int(sys.argv[3])
    b = int(sys.argv[4])
    c = float(sys.argv[5])
    out = sys.argv[6]
    virtual_center = (a, b, c)

    with open(data_dir + 'pdbs_train.txt') as file:
        pdbs_train = set(file.read().splitlines())

    # Find alignments between PDBs of the training set
    alignments = []
    with open(pairfile) as file:
        for line in file:
            sid1, sid2, cigar_string = line.rstrip('\n').split()

            if sid1 in pdbs_train and sid2 in pdbs_train:
                #print(' '.join((sid1, sid2, cigar_string)))
                alignments.append((sid1, sid2, cigar_string))

    # Not needed execept to exactly reproduce result
    random.Random(123).shuffle(alignments)

    xy = []  # (n x 10, n x 10)
    for sid1, sid2, cigar_string in alignments:
        xy.append(align_features(pdb_dir, virtual_center, sid1, sid2, cigar_string))

    # Write features to disc
    x_feat = np.vstack([x for x, y in xy])
    y_feat = np.vstack([y for x, y in xy])
    idx = np.arange(len(x_feat))
    np.random.RandomState(123).shuffle(idx)

    np.save(out, np.dstack([x_feat[idx], y_feat[idx]]))

