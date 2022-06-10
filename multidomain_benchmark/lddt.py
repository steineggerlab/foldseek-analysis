# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Code was adjusted by Martin Steinegger to support alignment comparison


"""lDDT protein distance score."""
import argparse

import numpy as np
import re

def lddt(target_points,
         query_points,
         aligned_pairs,
         cutoff=15.,
         per_residue=False):
  """Measure (approximate) lDDT for a batch of coordinates.

  lDDT reference:
  Mariani, V., Biasini, M., Barbato, A. & Schwede, T. lDDT: A local
  superposition-free score for comparing protein structures and models using
  distance difference tests. Bioinformatics 29, 2722–2728 (2013).

  lDDT is a measure of the difference between the true distance matrix and the
  distance matrix of the predicted points.  The difference is computed only on
  points closer than cutoff *in the true structure*.

  This function does not compute the exact lDDT value that the original paper
  describes because it does not include terms for physical feasibility
  (e.g. bond length violations). Therefore this is only an approximate
  lDDT score.

  Args:
    target_points: (batch, length, 3) array of predicted 3D points
    target_points: (batch, length, 3) array of true 3D points
    true_points_mask: (batch, length, 1) binary-valued float array.  This mask
      should be 1 for points that exist in the true points.
    cutoff: Maximum distance for a pair of points to be included
    per_residue: If true, return score for each residue.  Note that the overall
      lDDT is not exactly the mean of the per_residue lDDT's because some
      residues have more contacts than others.

  Returns:
    An (approximate, see above) lDDT score in the range 0-1.
  """

  assert len(target_points.shape) == 3
  assert target_points.shape[-1] == 3
  assert len(query_points.shape) == 3

  # Compute true and predicted distance matrices.
  dmat_query = np.sqrt(np.sum(
      (query_points[:, :, None] - query_points[:, None, :])**2, axis=-1))
  dists_to_score = (
      (dmat_query < cutoff).astype(np.float32) *
      (1. - np.eye(dmat_query.shape[1]))  # Exclude self-interaction.
  )

  # compute aligned distance matrix
  aligned_query_points=query_points[0, aligned_pairs[:, :, 0]]
  aligned_target_points=target_points[0, aligned_pairs[:, :, 1]]
  dmat_aligned_target = np.sqrt(np.sum(
      (aligned_target_points[:, :, None] -
       aligned_target_points[:, None, :])**2, axis=-1))
  dmat_aligned_query = np.sqrt(np.sum(
      (aligned_query_points[:, :, None] -
       aligned_query_points[:, None, :])**2, axis=-1))
  aligned_dists_to_score = (
      (dmat_aligned_query < cutoff).astype(np.float32) *
      (1. - np.eye(dmat_aligned_query.shape[1]))  # Exclude self-interaction.
  )

  dist_l1 = np.abs(dmat_aligned_query - dmat_aligned_target)

  # True lDDT uses a number of fixed bins.
  # We ignore the physical plausibility correction to lDDT, though.
  score = 0.25 * ((dist_l1 < 0.5).astype(np.float32) +
                  (dist_l1 < 1.0).astype(np.float32) +
                  (dist_l1 < 2.0).astype(np.float32) +
                  (dist_l1 < 4.0).astype(np.float32))
  # set diagonal of score to 0 (no self-interaction)
  score = score * (1. - np.eye(score.shape[1]))

  # Normalize over the appropriate axes.
  reduce_axes = (-1,) if per_residue else (-2, -1)
  norm = 1. / (np.sum(dists_to_score, axis=reduce_axes))
  norm_aligned = norm[0, aligned_pairs[:, :, 0]]

  score = norm_aligned * (np.sum(score * aligned_dists_to_score, axis=reduce_axes))
  return score

### Adjusted DeepMind code is over here ###

# read in pdb file without using bio.pdb
def read_ca_from_pdb(pdb_file):
  """Reads a PDB file and returns a numpy array of CA coordinates."""
  with open(pdb_file) as f:
    lines = f.readlines()
  coords = []
  for line in lines:
    if line.startswith('ATOM'):
      if line[12:16].strip() == 'CA':
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
  return np.array([coords])

# convert alignment cigar string to list of tuples of aligned residues positions
def cigar_to_tuples(cigar, query_start, target_start):
  """Convert a CIGAR string to a list of tuples of aligned residues.

  Args:
    cigar: A CIGAR string.
    query_start: The start position of the query sequence.
    target_start: The start position of the target sequence.

  Returns:
    A list of tuples of aligned residues.
  """
  tuples = []
  query_pos = query_start
  target_pos = target_start
  for length, op in re.findall(r'(\d+)([MDI])', cigar):
    if op == 'M':
      for i in range(int(length)):
        tuples.append((query_pos, target_pos))
        query_pos += 1
        target_pos += 1
    elif op == 'D':
      target_pos += int(length)
    elif op == 'I':
      query_pos += int(length)
  return np.array([tuples], dtype=np.int32)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # parse three arguments alignment file and query and pdb folder
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment', type=str, default='', help='alignment file')
    parser.add_argument('--query', type=str, default='', help='query pdb file')
    parser.add_argument('--target', type=str, default='', help='target pdb file')
    args = parser.parse_args()
    cutoff = 15.0
    per_residue = True
    print(args.query)
    # iterative over tsv file with five columns: query_id, target_id, query_start, target_start, cigar
    with open(args.alignment, 'r') as f:
        for line in f:
            query_id, target_id, query_start, target_start, cigar = line.split('\t')
            query_start = int(query_start) - 1
            target_start = int(target_start) - 1
            cigar = cigar.strip()
            query_pos = read_ca_from_pdb(f"{args.query}/{query_id}")
            target_pos = read_ca_from_pdb(f"{args.target}/{target_id}")
            aligned = cigar_to_tuples(cigar, query_start, target_start)
            per_residue_plddt = lddt(target_pos, query_pos, aligned, cutoff, per_residue)
            # join alignment and plddt scores
            query_lddt_info = ""
            for i,tuple in enumerate(aligned[0]):
                query_lddt_info += f"{tuple[0]}:{per_residue_plddt[0][i]:.3f},"
            print(f"{query_id}\t{target_id}\t{query_lddt_info[:-1]}")
