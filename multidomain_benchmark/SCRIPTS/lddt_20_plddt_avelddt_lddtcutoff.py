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
import os

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
  plddt = []
  for line in lines:
    if line.startswith('ATOM'):
      if line[12:16].strip() == 'CA':
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        plddt.append(float(line[61:66]))
  return np.array([coords]), plddt

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

# if lddt > threshold, increase TP counter, of <= threshold increase FP counter for that column
# only evaluate positions where query has plddt >= 70
def count_TP_per_column(query_count, plddt, lddt_per_coord, ave, lddt_cutoff):
     for lddt in lddt_per_coord:
         if plddt[lddt[0]] >= 70 and ave <= lddt_cutoff:
             if ave >= args.tp_threshold and query_count[lddt[0]][1] < 1:
                 query_count[lddt[0]][0] = query_count[lddt[0]][0] + 1
             elif ave <= args.fp_threshold:
                 query_count[lddt[0]][1]  = query_count[lddt[0]][1] + 1

def count_stuff_per_query(all_counts, normaliser):
    list = []
    for bla in all_counts:
        list.append(bla[0])
    n = 1
    per_query = ''
    while n <= max(list):
        per_query += str(n)
        per_query += ','
        per_query += str(sum(i >= n for i in list)/normaliser)
        n += 1
    print(per_query)

def count_stuff(all_counts, normaliser):
    list = []
    for count in all_counts:
        for bla in count:
            list.append(bla[0])
    n = 1
    while n <= max(list):
        print(n,sum(i >= n for i in list)/normaliser)
        n += 1

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # parse three arguments alignment file and query and pdb folder
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment', type=str, default='', help='alignment file')
    parser.add_argument('--query', type=str, default='', help='query pdb file')
    parser.add_argument('--target', type=str, default='', help='target pdb file')
    parser.add_argument('--tp-threshold', type=float, default=0.8, help='true positive threshold')
    parser.add_argument('--fp-threshold', type=float, default=0.2, help='false positive threshold')
    parser.add_argument('--lddt-cutoff', type=float, default=0.7)
    args = parser.parse_args()
    cutoff = 15.0
    per_residue = True
    all_counts = []
    query  = ''
    result = ''
    hit_counter = 0
    query_total = 0
    query_total_test = 0
    # get full query length
    with os.scandir(args.query) as it:
        for entry in it:
            query_pos, plddt = read_ca_from_pdb(entry)
            query_len_plddt = sum(i >= 70 for i in plddt)
            query_len = len(query_pos[0])
            query_total = query_total + query_len_plddt

    # iterative over tsv file with five columns: query_id, target_id, query_start, target_start, cigar
    with open(args.alignment, 'r') as f:
        for line in f:
            lddt_per_coord = []
            average_lddt = 0
            query_id, target_id, query_start, target_start, cigar = line.split('\t')
            query_start = int(query_start) - 1
            target_start = int(target_start) - 1
            cigar = cigar.strip()
            if query != query_id:
                query_pos, plddt_query = read_ca_from_pdb(f"{args.query}/{query_id}")
                query_len = len(query_pos[0])
            hit_counter = hit_counter + 1
            if query != query_id and query != '':
                all_counts.append(query_count)
                query_lddt_count = ''
                hit_counter = 0
                for i,tuple in enumerate(query_count):
                    query_lddt_count += f"{i}:{query_count[i][0]},"
                    result += str(query_count[i][0]) + ','
            if query == '' or query != query_id:
                query_count = [[0 for col in range(2)] for row in range(query_len)]
                query = query_id
            if hit_counter > 19:
                continue
            else:
                target_pos, plddt_target = read_ca_from_pdb(f"{args.target}/{target_id}")
                aligned = cigar_to_tuples(cigar, query_start, target_start)
                per_residue_lddt = lddt(target_pos, query_pos, aligned, cutoff, per_residue)
                average_lddt = sum(per_residue_lddt[0])/len(per_residue_lddt[0])
                # join alignment and plddt scores
                query_lddt_info = ""
                for i,tuple in enumerate(aligned[0]):
                    query_lddt_info += f"{tuple[0]}:{per_residue_lddt[0][i]:.3f},"
                    lddt_per_coord.append([tuple[0], per_residue_lddt[0][i]])
                count_TP_per_column(query_count, plddt_query, lddt_per_coord, average_lddt, args.lddt_cutoff)
    query_lddt_count = ''
    for i,tuple in enumerate(query_count):
        query_lddt_count += f"{i}:{query_count[i][0]},"
        result += str(query_count[i][0]) + ','
    all_counts.append(query_count)
    final_count = count_stuff(all_counts, query_total)
