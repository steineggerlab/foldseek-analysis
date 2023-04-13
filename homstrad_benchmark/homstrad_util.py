import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import warnings
import random
import pickle
import re
import os
import os.path as osp
import copy
import itertools


def pairs_from_aln(start1, aln1, start2, aln2):
    i1, i2 = start1, start2
    pairs = []
    for letter1, letter2 in zip(aln1, aln2):
        if letter1 not in ('-', '/') and letter2 not in ('-', '/'):
            pairs.append((i1, i2))
        if letter1 not in ('-', '/'):
            i1 += 1
        if letter2 not in ('-', '/'):
            i2 += 1
    return pairs

def parse_homstrad_ali(path):
    with open(path) as file:
        lines = file.readlines()
    aln = []
    
    while True:
        line = lines.pop(0)
        
        if line.startswith('C;'):
            pass
        elif not line.rstrip('\n'):
            pass
        elif line.startswith('>'):
            name = line[4:-1]
            aln.append((name, []))
            info_line = lines.pop(0)
            assert info_line.startswith('structure'), info_line + ' ' + path
        else:
            aln[-1][1].append(line)
            
        if not lines:
            break

    res = []
    for name, lst in aln:
        assert lst
        assert lst[-1].rstrip('\n').endswith('*'), f'{name} might be wrong'
        res.append((name, ''.join(line.rstrip('*\n') for line in lst)))
    return res

def parse_test_ali(path):
    aln = {}
    for name1, name2, start1, start2, qaln, taln in (i.split()[:6] for i in open(path)):
        if name1 != name2:
            aln[(name1, name2)] = (int(start1), qaln, int(start2), taln)
    return aln

def check_alignment(fam_name_fn, family, allvsall=False):
    aln = parse_test_ali(fam_name_fn(family))
    ref_aln = parse_homstrad_ali(f'homstrad_db/{family}/{family}.ali')

    res = []
    log = ''
    
    name1, ref_aln1 = ref_aln[0]
    name2, ref_aln2 = ref_aln[-1]

    if (name1, name2) not in aln:
        res.append([name1, name2, None, None, family, None])
        return res
    start1, aln1, start2, aln2 = predicted_aln = aln[(name1, name2)]

    ref_pairs = pairs_from_aln(0, ref_aln1, 0, ref_aln2)  # here
    pairs = pairs_from_aln(start1, aln1, start2, aln2)

    sensitivity = len([p for p in ref_pairs if p in pairs]) / len(ref_pairs)
    accuracy = len([p for p in pairs if p in ref_pairs]) / len(pairs)

    res.append([name1, name2, sensitivity, accuracy, family, predicted_aln])

    # Logging
    correct_pairs = [p for p in pairs if p in ref_pairs]
    correct_pos1 = {p[0] for p in correct_pairs}
    correct_pos2 = {p[1] for p in correct_pairs}
    refaln_pos1 = [p[0] for p in ref_pairs]
    refaln_pos2 = [p[1] for p in ref_pairs]
    seq1, between, seq2 = '', '', ''
    assert len(aln1) == len(aln2)
    i = start1
    for letter in aln1:
        if letter not in ('-', '/'):
            if i in correct_pos1:
                seq1 += str(refaln_pos1.index(i) % 10)  # 'o'
            else:
                if i in refaln_pos1:
                    #seq1 += 'x'
                    seq1 += str(refaln_pos1.index(i) % 10)  # 'x'
                else:
                    seq1 += 'n'  # FP
            i += 1
        else:
            seq1 += '-'
    i = start2
    for letter in aln2:
        if letter not in ('-', '/'):
            if i in correct_pos2:
                #seq2 += 'o'
                seq2 += str(refaln_pos2.index(i) % 10)  # 'o'
                between += '|'
            else:                        
                between += ' '
                if i in refaln_pos2:
                    seq2 += str(refaln_pos2.index(i) % 10)  # 'x'
                else:
                    seq2 += 'n'  # FP
            i += 1
        else:
            seq2 += '-'
            between += ' '

    log += f'> {name1} - {name2} (sensitivity: {sensitivity:.2f}, precision: {accuracy:.2f})\n\n'
    for i in range((len(aln1) + 99) // 100):
        log += seq1[i*100:(i+1)*100] + '\n' + between[i*100:(i+1)*100] + '\n' + seq2[i*100:(i+1)*100] + '\n\n'

    log += '*' * 100 + '\n\n'

                
    with open(fam_name_fn(family) + '.log', 'w') as file:
        file.write(log)
            
    return res

def check_alignments(fam_name_fn):
    res = []
    for family in tqdm(os.listdir('homstrad_db/')):
        res += check_alignment(fam_name_fn, family)        
    return res

def stats(res):
    found_ration = len([r for r in res if r[2] is not None]) / len(res)
    sensitivity = np.mean([r[2] or 0 for r in res])
    precision = np.mean([r[3] for r in res if r[3] is not None])
    return found_ration, sensitivity, precision

def seq_id(aln_a, aln_b):
    cnt = 0
    for a,b in zip(aln_a, aln_b):
        if a == '-' and b == '-':
            continue
        else:
            if a == b:
                cnt += 1
    return cnt

def read_fasta(fn):
    name2ss = {}
    with open(fn) as file:
        for line in file:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:]
            else:
                name2ss[name] = line
    return name2ss

def cnt_identities_hmstrd(aa_seq_aln1, aa_seq_aln2, ss_seq1, ss_seq2):
    aa_seq1 = aa_seq_aln1.replace('-', '').replace('/', '')
    aa_seq2 = aa_seq_aln2.replace('-', '').replace('/', '')
    
    # AA matches
    cnt_aa = seq_id(aa_seq_aln1, aa_seq_aln2)

    # 3Di matches
    aln = pairs_from_aln(0, aa_seq_aln1, 0, aa_seq_aln2)
    assert len(aa_seq1) == len(ss_seq1), f'{len(aa_seq1)} != {len(ss_seq1)}'
    assert len(aa_seq2) == len(ss_seq2)
    cnt_ss = 0
    for a, b in aln:
        if ss_seq1[a] == ss_seq2[b]:
            cnt_ss += 1
                
    return cnt_aa, cnt_ss

def read_submat(fn):
    import numpy as np
    lines = [line for line in open(fn).read().splitlines() if not line.startswith('#')]
    alphabet = lines[0].split()
    sub_mat = np.array([[float(i) for i in line[1:].split()] for line in lines[1:]])
    return sub_mat, alphabet



def rmsd(name1, name2, seq_a, seq_b, start1=1, start2=1, glob=True):
    if '../training' not in sys.path:  # assumption about cwd!
        sys.path.append('../training')
    import extract_pdb_features
    coords1, _ = extract_pdb_features.get_coords_from_pdb(name1)
    coords2, _ = extract_pdb_features.get_coords_from_pdb(name2)

    if glob:
        _seq1_len = len(seq_a.replace('-', '').replace('/', ''))
        assert _seq1_len == coords1.shape[0], f'{_seq1_len} vs {coords1.shape[0]}'
        assert len(seq_b.replace('-', '').replace('/', '')) == coords2.shape[0]

    aln_a, aln_b = [], []
    cnt_a, cnt_b, naln = start1-1, start2-1, 0
    for a,b in zip(seq_a, seq_b):
        if a != '-' or b != '-':
            if a == '-':
                cnt_b += 1
            elif b == '-':
                cnt_a += 1
            else:  # match
                # _d = np.sqrt(np.square(coords1[cnt_a, 0:3] - coords2[cnt_b, 0:3]).sum())
                # print(f'{cnt_a}: {a} | {cnt_b}: {b} => {_d:.2f} {coords1[cnt_a, 0:3]}, {coords2[cnt_b, 0:3]}')
                aln_a.append(coords1[cnt_a, 0:3])
                aln_b.append(coords2[cnt_b, 0:3])
                cnt_a += 1
                cnt_b += 1
                naln += 1

    rmsd = np.sqrt(np.square(np.linalg.norm(np.array(aln_a) - np.array(aln_b), axis=1)).mean())
    return rmsd, naln

def format_aln(q_name, t_name, q_aln, t_aln, q_start, t_start,
               line_len=90, include_middle=True, include_name=True, draw_plus=lambda a,b: False):
    assert len(q_aln) == len(t_aln)
    import math
    s = ''
    q_cnt, t_cnt = q_start, t_start
    for i in range(math.ceil(len(q_aln)/line_len)):
        q_line = q_aln[i*line_len:(i+1)*line_len]
        t_line = t_aln[i*line_len:(i+1)*line_len]
        diff = ''
        for a, b in zip(q_line, t_line):
            if a == b and a != '-':
                diff += '|'
            elif '-' not in (a, b) and draw_plus(a,b):
                diff += '+'
            else:
                diff += ' '
        if include_name:
            s += f'{q_name:<8}' + f'{q_cnt:>3}' + ' '
        s += q_line + '\n'
        if include_middle:
            if include_name:
                s += ' ' * 12
            s += diff + '\n'
        if include_name:
            s += f'{t_name:<8}' + f'{t_cnt:>3}' + ' '
        s += t_line + '\n'
        s += '\n'
        q_cnt += len(q_line.replace('-', ''))
        t_cnt += len(t_line.replace('-', ''))
    return s

def exchange_seq(aln_str, seq, require_same_len=True):
    if require_same_len:
        aln_res_no = len(aln_str.replace('-', ''))
        assert aln_res_no == len(seq), f'{aln_res_no} != {len(seq)}'
    s, cnt = '', 0
    for x in aln_str:
        if x == '-':
            s += '-'
        else:
            s += seq[cnt]
            cnt += 1
    return s