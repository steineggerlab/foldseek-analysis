from time import time
import argparse
from numpy import ndarray

import pickle
import typing
from dataclasses import dataclass
from pathlib import Path

import numba as nb
import numpy as np
import prody as pd
import typer
from geometricus import Structure, MomentInvariants, SplitType, GeometricusEmbedding
from scipy.spatial.distance import pdist, squareform
from copy import deepcopy

from caretta import (
    dynamic_time_warping as dtw,
    neighbor_joining as nj,
    score_functions,
    superposition_functions,
    feature_extraction,
    helper,
)


def alignment_to_numpy(alignment):
    aln_np = {}
    for n in alignment:
        aln_seq = []
        index = 0
        for a in alignment[n]:
            if a == "-":
                aln_seq.append(-1)
            else:
                aln_seq.append(index)
                index += 1
        aln_np[n] = np.array(aln_seq)
    return aln_np


@nb.njit
def tm_score(coords_1, coords_2, l1, l2):
    d1 = 1.24 * (l1 - 15) ** 1 / 3 - 1.8
    d2 = 1.24 * (l2 - 15) ** 1 / 3 - 1.8
    sum_1 = 0
    sum_2 = 0
    for i in range(coords_1.shape[0]):
        sum_1 += 1 / (1 + (np.sum(coords_1[i] - coords_2[i]) / d1)**2)
        sum_2 += 1 / (1 + (np.sum(coords_1[i] - coords_2[i]) / d2)**2)
    t1 = (1 / l1) * sum_1
    t2 = (1 / l2) * sum_2
    return max(t1, t2)


@nb.njit
def get_common_coordinates(
    coords_1: np.ndarray, coords_2: np.ndarray, aln_1: np.ndarray, aln_2: np.ndarray
) -> typing.Tuple[np.ndarray, np.ndarray]:
    """
    Return coordinate positions aligned in both coords_1 and coords_2
    """
    assert aln_1.shape == aln_2.shape
    pos_1, pos_2 = helper.get_common_positions(aln_1, aln_2)
    return coords_1[pos_1], coords_2[pos_2]


@nb.njit
def get_mean_coords(
    aln_1, coords_1: np.ndarray, aln_2, coords_2: np.ndarray
) -> np.ndarray:
    """
    Mean of two coordinate sets (of the same shape)
    """
    mean_coords = np.zeros((aln_1.shape[0], coords_1.shape[1]))
    for i, (x, y) in enumerate(zip(aln_1, aln_2)):
        if x == -1:
            mean_coords[i] = coords_2[y]
        elif y == -1:
            mean_coords[i] = coords_1[x]
        else:
            mean_coords[i] = np.array(
                [
                    np.nanmean(np.array([coords_1[x, d], coords_2[y, d]]))
                    for d in range(coords_1.shape[1])
                ]
            )
    return mean_coords


def get_mean_weights(
    weights_1: np.ndarray, weights_2: np.ndarray, aln_1: np.ndarray, aln_2: np.ndarray
) -> np.ndarray:
    mean_weights = np.zeros(aln_1.shape[0])
    for i, (x, y) in enumerate(zip(aln_1, aln_2)):
        if not x == -1:
            mean_weights[i] += weights_1[x]
        if not y == -1:
            mean_weights[i] += weights_2[y]
    return mean_weights


@nb.njit
def get_pairwise_alignment(
    coords_1,
    coords_2,
    gamma,
    gap_open_penalty: float,
    gap_extend_penalty: float,
    weights_1: np.ndarray,
    weights_2: np.ndarray,
    n_iter=3,
):
    score_matrix = score_functions.make_score_matrix(
        np.hstack((coords_1, weights_1)),
        np.hstack((coords_2, weights_2)),
        score_functions.get_caretta_score,
        gamma,
    )
    dtw_aln_array_1, dtw_aln_array_2, dtw_score = dtw.dtw_align(
        score_matrix, gap_open_penalty, gap_extend_penalty
    )
    for i in range(n_iter):
        pos_1, pos_2 = helper.get_common_positions(dtw_aln_array_1, dtw_aln_array_2)
        common_coords_1, common_coords_2 = coords_1[pos_1], coords_2[pos_2]
        (
            c1,
            c2,
            common_coords_2,
        ) = superposition_functions.paired_svd_superpose_with_subset(
            coords_1, coords_2, common_coords_1, common_coords_2
        )
        score_matrix = score_functions.make_score_matrix(
            np.hstack((c1, weights_1)),
            np.hstack((c2, weights_2)),
            score_functions.get_caretta_score,
            gamma,
        )
        aln_1, aln_2, score = dtw.dtw_align(
            score_matrix, gap_open_penalty, gap_extend_penalty
        )
        if score > dtw_score:
            coords_1 = c1
            coords_2 = c2
            dtw_score = score
            dtw_aln_array_1 = aln_1
            dtw_aln_array_2 = aln_2
        else:
            break
    return dtw_aln_array_1, dtw_aln_array_2, dtw_score, coords_1, coords_2


@dataclass
class OutputFiles:
    fasta_file: Path = Path("./result.fasta")
    pdb_folder: Path = Path("./result_pdb/")
    cleaned_pdb_folder: Path = Path("./cleaned_pdb")
    feature_file: Path = Path("./result_features.pkl")
    class_file: Path = Path("./result_class.pkl")


DEFAULT_SUPERPOSITION_PARAMETERS = {
    # must-have
    "gap_open_penalty": 0.0,
    "gap_extend_penalty": 0.0,
    "gamma": 0.03,
    # changes per superposition_function
    "split_type": "KMER",
    "split_size": 20,
    "scale": True,
    "gamma_moment": 0.6,
    "n_iter": 3,
}


@dataclass
class StructureMultiple:
    """
    Class for multiple structure alignment

    Constructor Arguments
    ---------------------
    structures
        list of protein_utility.Structure objects
    superposition_parameters
        dictionary of parameters to pass to the superposition function
    superposition_function
        a function that takes two coordinate sets as input and superposes them
        returns a score, superposed_coords_1, superposed_coords_2
    score_function
        a function that takes two paired coordinate sets (same shape) and returns a score
    consensus_weight
        weights the effect of well-aligned columns on the progressive alignment
    final_structures
        makes the progressive alignment tree of structures, last one is the consensus structure of the multiple alignment
    tree
        neighbor joining tree of indices controlling alignment - indices beyond len(structures) refer to intermediate nodes
    branch_lengths
    alignment
        indices of aligning residues from each structure, gaps are -1s
    """

    structures: typing.List[Structure]
    sequences: typing.Dict[str, str]
    superposition_parameters: typing.Dict[str, typing.Any]
    superposition_function: typing.Callable[
        [
            np.ndarray,
            np.ndarray,
            dict,
            typing.Callable[[np.ndarray, np.ndarray], float],
        ],
        typing.Tuple[float, np.ndarray, np.ndarray],
    ] = lambda x, y: (0, x, y)
    mean_function: typing.Callable[
        [np.ndarray, np.ndarray, np.ndarray, np.ndarray], np.ndarray
    ] = get_mean_coords
    consensus_weight: float = 1.0
    pairwise_distance_matrix: typing.Union[None, np.ndarray] = None
    reference_structure_index: typing.Union[None, int] = None
    final_structures: typing.Union[None, typing.List[Structure]] = None
    final_consensus_weights: typing.Union[None, typing.List[np.ndarray]] = None
    tree: typing.Union[None, np.ndarray] = None
    branch_lengths: typing.Union[None, np.ndarray] = None
    alignment: typing.Union[dict, None] = None
    output_folder: Path = Path("./caretta_results")
    features: typing.Union[dict, None] = None

    @staticmethod
    def align_from_pdb_files(
        input_pdb: typing.Union[typing.List[str], Path, str],
        gap_open_penalty: float = 1.0,
        gap_extend_penalty: float = 0.01,
        consensus_weight: float = 1.0,
        full: bool = False,
        output_folder: typing.Union[str, Path] = Path("../caretta_results"),
        num_threads: int = 20,
        write_fasta: bool = False,
        write_pdb: bool = False,
        write_features: bool = False,
        only_dssp: bool = True,
        write_class: bool = False,
        write_matrix: bool = False,
        verbose: bool = True,
    ):
        """
        Caretta aligns protein structures and can output a sequence alignment, superposed PDB files,
        a set of aligned feature matrices and a class with intermediate structures made during progressive alignment.

        Parameters
        ----------
        input_pdb
            Can be \n
            A list of PDB files,
            A list of PDB IDs,
            A folder with input protein files,
            A file which lists PDB filenames on each line,
            A file which lists PDB IDs on each line,
        gap_open_penalty
            default 1
        gap_extend_penalty
            default 0.01
        consensus_weight
            default 1
        full
            True =>  Uses all-vs-all pairwise Caretta alignment to make the distance matrix (much slower)
        output_folder
            default "caretta_results"
        num_threads
            Number of threads to use for feature extraction
        write_fasta
            True => writes alignment as fasta file (default True)
            writes to output_folder / result.fasta
        write_pdb
            True => writes all protein PDB files superposed by alignment (default True)
             writes to output_folder / superposed_pdb
        write_features
            True => extracts and writes aligned features as a dictionary of numpy arrays into a pickle file (default True)
            writes to output_folder / result_features.pkl
        only_dssp
            True => only DSSP features extracted
            False => all features
        write_class
            True => writes StructureMultiple class with intermediate structures and tree to pickle file (default True)
            writes to output_folder / result_class.pkl
        write_matrix
            True => writes distance matrix to text file (default False)
            writes to output_folder / matrix.mat
        verbose
            controls verbosity

        Returns
        -------
        StructureMultiple class
        """
        msa_class = StructureMultiple.from_pdb_files(
            input_pdb,
            superposition_parameters=DEFAULT_SUPERPOSITION_PARAMETERS,
            superposition_function=superposition_functions.moment_svd_superpose_function,
            consensus_weight=consensus_weight,
            output_folder=output_folder,
            verbose=verbose,
        )
        if len(msa_class.structures) > 2:
            if full:
                msa_class.make_pairwise_dtw_matrix(
                    gap_open_penalty, gap_extend_penalty, verbose=verbose
                )
            else:
                msa_class.make_pairwise_shape_matrix(
                    verbose=verbose
                )
            msa_class.align(gap_open_penalty, gap_extend_penalty, verbose=verbose)
        else:
            msa_class.reference_structure_index = 0
            msa_class.align(
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                verbose=verbose,
            )

        msa_class.write_files(
            write_fasta=write_fasta,
            write_pdb=write_pdb,
            write_features=write_features,
            write_class=write_class,
            write_matrix=write_matrix,
            num_threads=num_threads,
            only_dssp=only_dssp,
            verbose=verbose,
        )
        return msa_class

    @classmethod
    def from_pdb_files(
        cls,
        input_pdb,
        superposition_parameters=DEFAULT_SUPERPOSITION_PARAMETERS,
        superposition_function=superposition_functions.moment_svd_superpose_function,
        consensus_weight=1.0,
        output_folder=Path("./caretta_results"),
        verbose: bool = False,
    ):
        """
        Makes a StructureMultiple object from a list of pdb files/names or a folder of pdb files

        Parameters
        ----------
        input_pdb
            list of pdb files/names or a folder containing pdb files
        superposition_parameters
            parameters to give to the superposition function
        superposition_function
            a function that takes two coordinate sets as input and superposes them
            returns a score, superposed_coords_1, superposed_coords_2
        consensus_weight
            weights the effect of well-aligned columns on the progressive alignment
        output_folder
        verbose

        Returns
        -------
        StructureMultiple object (unaligned)
        """
        output_folder = Path(output_folder)
        if not output_folder.exists():
            output_folder.mkdir()

        cleaned_pdb_folder = output_folder / "cleaned_pdb"
        if not cleaned_pdb_folder.exists():
            cleaned_pdb_folder.mkdir()
        pdb_files = helper.parse_pdb_files_and_clean(input_pdb, cleaned_pdb_folder)
        if verbose:
            typer.echo(f"Found {len(pdb_files)} PDB files")

        structures = []
        sequences = {}
        for pdb_file in pdb_files:
            pdb_name = Path(pdb_file).stem
            protein = pd.parsePDB(str(pdb_file)).select("protein and calpha")
            coordinates = protein.getCoords()
            structures.append(Structure(pdb_name, coordinates.shape[0], coordinates))
            sequences[pdb_name] = protein.getSequence()
        msa_class = StructureMultiple(
            structures,
            sequences,
            superposition_parameters,
            superposition_function,
            consensus_weight=consensus_weight,
            output_folder=output_folder,
        )
        return msa_class

    @classmethod
    def from_coordinates(
        cls,
        names: typing.List[str],
        coordinates_list: typing.List[np.ndarray],
        sequences: typing.List[str],
        superposition_parameters: dict,
        superposition_function=superposition_functions.moment_multiple_svd_superpose_function,
        consensus_weight=1.0,
        output_folder=Path("./caretta_results"),
    ):
        """
        Makes a StructureMultiple object from a list of coordinates

        Parameters
        ----------
        names
        coordinates_list
        sequences
        superposition_parameters
            parameters to give to the superposition function
        superposition_function
            a function that takes two coordinate sets as input and superposes them
            returns a score, superposed_coords_1, superposed_coords_2
        consensus_weight
            weights the effect of well-aligned columns on the progressive alignment
        output_folder

        Returns
        -------
        StructureMultiple object (unaligned)
        """
        output_folder = Path(output_folder)
        if not output_folder.exists():
            output_folder.mkdir()
        sequences = {n: s for n, s in zip(names, sequences)}
        structures = []
        for name, coordinates in zip(names, coordinates_list):
            structures.append(Structure(name, coordinates.shape[0], coordinates))
        msa_class = StructureMultiple(
            structures,
            sequences,
            superposition_parameters,
            superposition_function,
            consensus_weight=consensus_weight,
            output_folder=output_folder,
        )
        return msa_class

    def get_pairwise_alignment(
        self,
        coords_1,
        coords_2,
        gap_open_penalty: float,
        gap_extend_penalty: float,
        weight=False,
        weights_1=None,
        weights_2=None,
        n_iter=3,
        verbose: bool = False,
    ):
        """
        Aligns coords_1 to coords_2 by first superposing and then running dtw on the score matrix of the superposed coordinate sets
        Parameters
        ----------
        coords_1
        coords_2
        gap_open_penalty
        gap_extend_penalty
        weight
        weights_1
        weights_2
        n_iter
        verbose

        Returns
        -------
        alignment_1, alignment_2, score, superposed_coords_1, superposed_coords_2
        """
        _, coords_1, coords_2 = self.superposition_function(
            coords_1, coords_2, self.superposition_parameters
        )

        if weight:
            assert weights_1 is not None
            assert weights_2 is not None
            weights_1 = weights_1.reshape(-1, 1)
            weights_2 = weights_2.reshape(-1, 1)
        else:
            weights_1 = np.zeros((coords_1.shape[0], 1))
            weights_2 = np.zeros((coords_2.shape[0], 1))
        return get_pairwise_alignment(
            coords_1,
            coords_2,
            self.superposition_parameters["gamma"],
            gap_open_penalty,
            gap_extend_penalty,
            weights_1,
            weights_2,
            n_iter,
        )

    def make_pairwise_shape_matrix(
        self,
        resolution: typing.Union[float, np.ndarray] = 2.0,
        parameters: dict = None,
        metric="braycurtis",
        verbose: bool = False,
    ):
        """
        Makes an all vs. all matrix of distance scores between all the structures.

        Parameters
        ----------
        resolution
        parameters
            to use for making invariants
            needs to have
            num_split_types
            split_type_i
            split_size_i
        metric
            distance metric (accepts any metric supported by scipy.spatial.distance
        verbose
        Returns
        -------
        [n x n] distance matrix
        """
        if verbose:
            typer.echo("Calculating pairwise distances...")
        if parameters is None:
            parameters = dict(num_split_types=1, split_type_0="KMER", split_size_0=20)
        embedders = []
        for i in range(parameters["num_split_types"]):
            invariants = (
                MomentInvariants.from_coordinates(
                    s.name,
                    s.coordinates,
                    None,
                    split_size=parameters[f"split_size_{i}"],
                    split_type=SplitType[parameters[f"split_type_{i}"]],
                )
                for s in self.structures
            )
            embedders.append(
                GeometricusEmbedding.from_invariants(
                    invariants,
                    resolution=resolution,
                    protein_keys=[s.name for s in self.structures],
                )
            )
        distance_matrix = squareform(
            pdist(
                np.hstack([embedder.embedding for embedder in embedders]),
                metric=metric,
            )
        )
        self.reference_structure_index = np.argmin(
            np.median(distance_matrix, axis=0)
        )
        self.pairwise_distance_matrix = distance_matrix

    def make_pairwise_dtw_matrix(
        self,
        gap_open_penalty: float,
        gap_extend_penalty: float,
        invert=True,
        verbose: bool = False,
    ):
        """
        Makes an all vs. all matrix of distance (or similarity) scores between all the structures using pairwise alignment.

        Parameters
        ----------
        gap_open_penalty
        gap_extend_penalty
        invert
            if True returns distance matrix
            if False returns similarity matrix
        verbose
        Returns
        -------
        [n x n] matrix
        """
        if verbose:
            typer.echo("Calculating pairwise distances...")
        pairwise_matrix = np.zeros((len(self.structures), len(self.structures)))
        for i in range(pairwise_matrix.shape[0] - 1):
            for j in range(i + 1, pairwise_matrix.shape[1]):
                coords_1, coords_2 = (
                    self.structures[i].coordinates,
                    self.structures[j].coordinates,
                )
                (
                    dtw_aln_1,
                    dtw_aln_2,
                    score,
                    coords_1,
                    coords_2,
                ) = self.get_pairwise_alignment(
                    coords_1,
                    coords_2,
                    gap_open_penalty=gap_open_penalty,
                    gap_extend_penalty=gap_extend_penalty,
                    weight=False,
                )
                common_coords_1, common_coords_2 = get_common_coordinates(
                    coords_1, coords_2, dtw_aln_1, dtw_aln_2
                )
                pairwise_matrix[i, j] = score_functions.get_total_score(
                    common_coords_1,
                    common_coords_2,
                    score_functions.get_caretta_score,
                    self.superposition_parameters["gamma"],
                    True,
                )
                if invert:
                    pairwise_matrix[i, j] *= -1
        pairwise_matrix += pairwise_matrix.T
        self.reference_structure_index = np.argmin(
            np.median(pairwise_matrix, axis=0)
        )
        self.pairwise_distance_matrix = pairwise_matrix

    def align(
        self,
        gap_open_penalty,
        gap_extend_penalty,
        return_sequence: bool = True,
        verbose: bool = False,
    ) -> dict:
        """
        Makes a multiple structure alignment

        Parameters
        ----------
        pw_matrix
            pairwise similarity matrix to base the neighbor joining tree on
        gap_open_penalty
        gap_extend_penalty
        return_sequence
            if True returns sequence alignment
            else indices of aligning residues with gaps as -1s
        verbose
        Returns
        -------
        alignment = {name: indices of aligning residues with gaps as -1s}
        """
        if len(self.structures) == 2:
            coords_1, coords_2 = (
                self.structures[0].coordinates,
                self.structures[1].coordinates,
            )
            dtw_1, dtw_2, _, _, _ = self.get_pairwise_alignment(
                coords_1,
                coords_2,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                weight=False,
                n_iter=self.superposition_parameters["n_iter"],
                verbose=verbose,
            )
            self.alignment = {
                self.structures[0].name: dtw_1,
                self.structures[1].name: dtw_2,
            }
            if return_sequence:
                return self.make_sequence_alignment()
            else:
                return self.alignment
        assert self.pairwise_distance_matrix is not None
        assert self.pairwise_distance_matrix.shape[0] == len(self.structures)
        if verbose:
            typer.echo("Constructing neighbor joining tree...")
        tree, branch_lengths = nj.neighbor_joining(self.pairwise_distance_matrix)
        self.tree = tree
        self.branch_lengths = branch_lengths
        self.final_structures = [s for s in self.structures]
        self.final_consensus_weights = [
            np.full(
                (self.structures[i].coordinates.shape[0], 1),
                self.consensus_weight,
                dtype=np.float64,
            )
            for i in range(len(self.structures))
        ]
        msa_alignments = {
            s.name: {s.name: np.arange(s.length)} for s in self.final_structures
        }

        def make_intermediate_node(n1, n2, n_int):
            name_1, name_2 = (
                self.final_structures[n1].name,
                self.final_structures[n2].name,
            )

            name_int = f"int-{n_int}"
            n1_coords = self.final_structures[n1].coordinates
            n1_weights = self.final_consensus_weights[n1]
            n1_weights *= len(msa_alignments[name_2])
            n1_weights /= 2 * (
                len(msa_alignments[name_2]) + len(msa_alignments[name_1])
            )
            n2_coords = self.final_structures[n2].coordinates
            n2_weights = self.final_consensus_weights[n2]
            n2_weights *= len(msa_alignments[name_1])
            n2_weights /= 2 * (
                len(msa_alignments[name_2]) + len(msa_alignments[name_1])
            )
            (
                dtw_aln_1,
                dtw_aln_2,
                score,
                n1_coords,
                n2_coords,
            ) = self.get_pairwise_alignment(
                n1_coords,
                n2_coords,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                weight=True,
                weights_1=n1_weights,
                weights_2=n2_weights,
                n_iter=self.superposition_parameters["n_iter"],
            )
            n1_weights *= (
                2
                * (len(msa_alignments[name_2]) + len(msa_alignments[name_1]))
                / len(msa_alignments[name_2])
            )
            n2_weights *= (
                2
                * (len(msa_alignments[name_2]) + len(msa_alignments[name_1]))
                / len(msa_alignments[name_1])
            )
            msa_alignments[name_1] = {
                name: np.array([sequence[i] if i != -1 else -1 for i in dtw_aln_1])
                for name, sequence in msa_alignments[name_1].items()
            }
            msa_alignments[name_2] = {
                name: np.array([sequence[i] if i != -1 else -1 for i in dtw_aln_2])
                for name, sequence in msa_alignments[name_2].items()
            }
            msa_alignments[name_int] = {
                **msa_alignments[name_1],
                **msa_alignments[name_2],
            }

            mean_coords = self.mean_function(dtw_aln_1, n1_coords, dtw_aln_2, n2_coords)
            mean_weights = get_mean_weights(
                n1_weights, n2_weights, dtw_aln_1, dtw_aln_2
            )
            self.final_structures.append(
                Structure(name_int, mean_coords.shape[0], mean_coords)
            )
            self.final_consensus_weights.append(mean_weights)

        if verbose:
            with typer.progressbar(
                range(0, self.tree.shape[0] - 1, 2), label="Aligning"
            ) as progress:
                for x in progress:
                    node_1, node_2, node_int = (
                        self.tree[x, 0],
                        self.tree[x + 1, 0],
                        self.tree[x, 1],
                    )
                    assert self.tree[x + 1, 1] == node_int
                    make_intermediate_node(node_1, node_2, node_int)
        else:
            for x in range(0, self.tree.shape[0] - 1, 2):
                node_1, node_2, node_int = (
                    self.tree[x, 0],
                    self.tree[x + 1, 0],
                    self.tree[x, 1],
                )
                assert self.tree[x + 1, 1] == node_int
                make_intermediate_node(node_1, node_2, node_int)

        node_1, node_2 = self.tree[-1, 0], self.tree[-1, 1]
        make_intermediate_node(node_1, node_2, "final")
        alignment = {
            **msa_alignments[self.final_structures[node_1].name],
            **msa_alignments[self.final_structures[node_2].name],
        }
        self.alignment = alignment
        if return_sequence:
            return self.make_sequence_alignment()
        else:
            return alignment

    def make_sequence_alignment(self, alignment=None):
        sequence_alignment = {}
        if alignment is None:
            alignment = self.alignment
        for s in self.structures:
            sequence_alignment[s.name] = "".join(
                self.sequences[s.name][i] if i != -1 else "-" for i in alignment[s.name]
            )
        return sequence_alignment

    def write_files(
        self,
        write_fasta,
        write_pdb,
        write_features,
        write_class,
        write_matrix,
        only_dssp=True,
        num_threads=4,
        verbose: bool = False,
    ):
        if verbose and any(
            (write_fasta, write_pdb, write_pdb, write_class, write_matrix)
        ):
            typer.echo("Writing files...")
        if write_fasta:
            fasta_file = self.output_folder / "result.fasta"
            self.write_alignment(fasta_file)
            if verbose:
                typer.echo(
                    f"FASTA file: {typer.style(str(fasta_file), fg=typer.colors.GREEN)}",
                )
        if write_pdb:
            pdb_folder = self.output_folder / "superposed_pdbs"
            if not pdb_folder.exists():
                pdb_folder.mkdir()
            self.write_superposed_pdbs(pdb_folder, verbose=verbose)
            if verbose:
                typer.echo(
                    f"Superposed PDB files: {typer.style(str(pdb_folder), fg=typer.colors.GREEN)}"
                )
        if write_features:
            dssp_dir = self.output_folder / ".caretta_tmp"
            if not dssp_dir.exists():
                dssp_dir.mkdir()
            feature_file = self.output_folder / "result_features.pkl"
            names, self.features = self.get_aligned_features(
                str(dssp_dir), only_dssp=only_dssp, num_threads=num_threads
            )
            with open(feature_file, "wb") as f:
                pickle.dump((names, self.features), f)
            if verbose:
                typer.echo(
                    f"Aligned features: {typer.style(str(feature_file), fg=typer.colors.GREEN)}"
                )
        if write_class:
            class_file = self.output_folder / "result_class.pkl"
            with open(class_file, "wb") as f:
                pickle.dump(self, f)
            if verbose:
                typer.echo(
                    f"Class file: {typer.style(str(class_file), fg=typer.colors.GREEN)}"
                )
        if write_matrix:
            matrix_file = self.output_folder / "matrix.mat"
            helper.write_distance_matrix(
                [s.name for s in self.structures],
                self.pairwise_distance_matrix,
                matrix_file,
            )
            if verbose:
                typer.echo(
                    f"Distance matrix file: {typer.style(str(matrix_file), fg=typer.colors.GREEN)}"
                )

    def write_alignment(self, filename, alignments: dict = None):
        """
        Writes alignment to a fasta file
        """
        if alignments is None:
            alignments = self.alignment
        with open(filename, "w") as f:
            for key in alignments:
                sequence = "".join(
                    self.sequences[key][n] if n != -1 else "-" for n in alignments[key]
                )
                f.write(f">{key}\n{sequence}\n")

    def write_superposed_pdbs(
        self, output_pdb_folder, alignments: dict = None, verbose: bool = False
    ):
        """
        Superposes PDBs according to alignment and writes transformed PDBs to files
        (View with Pymol)

        Parameters
        ----------
        output_pdb_folder
        alignments
        verbose
        """
        if alignments is None:
            alignments = self.alignment
        output_pdb_folder = Path(output_pdb_folder)
        if not output_pdb_folder.exists():
            output_pdb_folder.mkdir()
        reference_name = self.structures[self.reference_structure_index].name
        core_indices = np.array(
            [
                i
                for i in range(len(alignments[reference_name]))
                if -1 not in [alignments[n][i] for n in alignments]
            ]
        )
        if verbose:
            typer.echo(
                f"{len(core_indices)} core positions in alignment of length {len(alignments[reference_name])}"
            )
        if len(core_indices) < len(alignments[reference_name]) // 2:
            if verbose:
                typer.echo(
                    typer.style(
                        "Core indices are < half of alignment length, superposing using reference structure "
                        "instead",
                        fg=typer.colors.RED,
                    )
                )
                typer.echo(
                    typer.style(
                        "Please inspect the distance matrix to split divergent protein groups",
                        fg=typer.colors.RED,
                    )
                )
            self.write_superposed_pdbs_reference(output_pdb_folder, alignments)
        else:
            self.write_superposed_pdbs_core(output_pdb_folder, alignments)

    def write_superposed_pdbs_core(self, output_pdb_folder, alignments):
        """
        Superposes PDBs according to core indices in alignment and writes transformed PDBs to files
        (View with Pymol)

        Parameters
        ----------
        alignments
        output_pdb_folder
        """
        reference_name = self.structures[self.reference_structure_index].name
        core_indices = np.array(
            [
                i
                for i in range(len(alignments[reference_name]))
                if -1 not in [alignments[n][i] for n in alignments]
            ]
        )
        reference_pdb = pd.parsePDB(
            str(
                self.output_folder
                / f"cleaned_pdb/{self.structures[self.reference_structure_index].name}.pdb"
            )
        )
        aln_ref = alignments[reference_name]
        ref_coords_core = (
            reference_pdb[helper.get_alpha_indices(reference_pdb)]
            .getCoords()
            .astype(np.float64)[np.array([aln_ref[c] for c in core_indices])]
        )
        ref_centroid = helper.nb_mean_axis_0(ref_coords_core)
        ref_coords_core -= ref_centroid
        transformation = pd.Transformation(np.eye(3), -ref_centroid)
        reference_pdb = pd.applyTransformation(transformation, reference_pdb)
        pd.writePDB(str(output_pdb_folder / f"{reference_name}.pdb"), reference_pdb)
        for i in range(1, len(self.structures)):
            name = self.structures[i].name
            pdb = pd.parsePDB(
                str(self.output_folder / f"cleaned_pdb/{self.structures[i].name}.pdb")
            )
            aln_name = alignments[name]
            common_coords_2 = (
                pdb[helper.get_alpha_indices(pdb)]
                .getCoords()
                .astype(np.float64)[np.array([aln_name[c] for c in core_indices])]
            )
            (
                rotation_matrix,
                translation_matrix,
            ) = superposition_functions.svd_superimpose(
                ref_coords_core, common_coords_2
            )
            transformation = pd.Transformation(rotation_matrix.T, translation_matrix)
            pdb = pd.applyTransformation(transformation, pdb)
            pd.writePDB(str(output_pdb_folder / f"{name}.pdb"), pdb)

    def write_superposed_pdbs_reference(self, output_pdb_folder, alignments):
        """
        Superposes PDBs according to reference structure and writes transformed PDBs to files
        (View with Pymol)

        Parameters
        ----------
        alignments
        output_pdb_folder
        """
        reference_name = self.structures[self.reference_structure_index].name
        reference_pdb = pd.parsePDB(
            str(
                self.output_folder
                / f"cleaned_pdb/{self.structures[self.reference_structure_index].name}.pdb"
            )
        )
        aln_ref = alignments[reference_name]
        reference_coords = (
            reference_pdb[helper.get_alpha_indices(reference_pdb)]
            .getCoords()
            .astype(np.float64)
        )
        pd.writePDB(str(output_pdb_folder / f"{reference_name}.pdb"), reference_pdb)
        for i in range(len(self.structures)):
            if i == self.reference_structure_index:
                continue
            name = self.structures[i].name
            pdb = pd.parsePDB(
                str(self.output_folder / f"cleaned_pdb/{self.structures[i].name}.pdb")
            )
            aln_name = alignments[name]
            common_coords_1, common_coords_2 = get_common_coordinates(
                reference_coords,
                pdb[helper.get_alpha_indices(pdb)].getCoords().astype(np.float64),
                aln_ref,
                aln_name,
            )
            (
                rotation_matrix,
                translation_matrix,
            ) = superposition_functions.svd_superimpose(
                common_coords_1, common_coords_2
            )
            transformation = pd.Transformation(rotation_matrix.T, translation_matrix)
            pdb = pd.applyTransformation(transformation, pdb)
            pd.writePDB(str(output_pdb_folder / f"{name}.pdb"), pdb)

    def get_aligned_features(
        self, dssp_dir, num_threads, alignment: dict = None, only_dssp: bool = True
    ) -> typing.Tuple[typing.List[str], typing.Dict[str, ndarray]]:
        """
        Get list of protein names and corresponding dict of aligned features
        """
        if alignment is None:
            alignment = self.make_sequence_alignment()

        pdb_files = [
            self.output_folder / "cleaned_pdb" / f"{s.name}.pdb"
            for s in self.structures
        ]
        features = feature_extraction.get_features_multiple(
            pdb_files,
            str(dssp_dir),
            num_threads=num_threads,
            only_dssp=only_dssp,
            force_overwrite=True,
        )
        feature_names = list(features[0].keys())
        aligned_features = {}
        alignment_length = len(alignment[self.structures[0].name])
        for feature_name in feature_names:
            if feature_name == "secondary":
                continue
            aligned_features[feature_name] = np.zeros(
                (len(self.structures), alignment_length)
            )
            aligned_features[feature_name][:] = np.nan
            for p in range(len(self.structures)):
                farray = features[p][feature_name]
                if "gnm" in feature_name or "anm" in feature_name:
                    farray = farray / np.nansum(farray ** 2) ** 0.5
                indices = [
                    i
                    for i in range(alignment_length)
                    if alignment[self.structures[p].name][i] != "-"
                ]
                aligned_features[feature_name][p, indices] = farray
        return [p.name for p in self.structures], aligned_features

    def superpose(self, alignments: dict = None):
        """
        Superposes structures to first structure using Kabsch superposition
        """
        if alignments is None:
            alignments = self.alignment
        reference_key = self.structures[self.reference_structure_index].name
        core_indices = np.array(
            [
                i
                for i in range(len(alignments[reference_key]))
                if "-" not in [alignments[n][i] for n in alignments]
            ]
        )
        if len(core_indices) < len(alignments[reference_key]) // 2:
            self.superpose_reference(alignments)
        else:
            self.superpose_core(alignments, core_indices)

    def superpose_core(self, alignments: dict = None, core_indices: np.ndarray = None):
        """
        Superposes structures to first structure according to core positions in alignment using Kabsch superposition
        """
        if alignments is None:
            alignments = self.alignment
        reference_key = self.structures[self.reference_structure_index].name
        if core_indices is None:
            core_indices = np.array(
                [
                    i
                    for i in range(len(alignments[reference_key]))
                    if "-" not in [alignments[n][i] for n in alignments]
                ]
            )
        aln_ref = alignments[reference_key]
        ref_coords = self.structures[self.reference_structure_index].coordinates[
            np.array([aln_ref[c] for c in core_indices])
        ]
        ref_centroid = helper.nb_mean_axis_0(ref_coords)
        ref_coords -= ref_centroid
        for i in range(len(self.structures)):
            if i == self.reference_structure_index:
                self.structures[i].coordinates -= ref_centroid
            else:
                aln_c = alignments[self.structures[i].name]
                common_coords_2 = self.structures[i].coordinates[
                    np.array([aln_c[c] for c in core_indices])
                ]
                (
                    rotation_matrix,
                    translation_matrix,
                ) = superposition_functions.paired_svd_superpose(
                    ref_coords, common_coords_2
                )
                self.structures[i].coordinates = superposition_functions.apply_rotran(
                    self.structures[i].coordinates, rotation_matrix, translation_matrix
                )

    def superpose_reference(self, alignments: dict = None):
        """
        Superposes structures to first structure according to reference structure using Kabsch superposition
        """
        if alignments is None:
            alignments = self.alignment
        reference_key = self.structures[self.reference_structure_index].name
        aln_ref = alignments[reference_key]
        for i in range(len(self.structures)):
            aln_c = alignments[self.structures[i].name]
            common_coords_1, common_coords_2 = get_common_coordinates(
                self.structures[self.reference_structure_index].coordinates,
                self.structures[i].coordinates,
                aln_ref,
                aln_c,
            )
            assert common_coords_1.shape[0] > 0
            (
                rotation_matrix,
                translation_matrix,
            ) = superposition_functions.paired_svd_superpose(
                common_coords_1, common_coords_2
            )
            self.structures[i].coordinates = superposition_functions.apply_rotran(
                self.structures[i].coordinates, rotation_matrix, translation_matrix
            )

    def make_rmsd_coverage_tm_matrix(
        self, alignments: dict = None, superpose_first: bool = True
    ):
        """
        Find RMSDs and coverages of the alignment of each pair of sequences

        Parameters
        ----------
        alignments
            if None uses self.alignment
        superpose_first
            if True then superposes all structures to first structure first

        Returns
        -------
        RMSD matrix, coverage matrix
        """
        if alignments is None:
            alignments = self.alignment
        num = len(self.structures)
        pairwise_rmsd_matrix = np.zeros((num, num))
        pairwise_rmsd_matrix[:] = np.nan
        pairwise_coverage = np.zeros((num, num))
        pairwise_coverage[:] = np.nan
        pairwise_tm = np.zeros((num, num))
        pairwise_tm[:] = np.nan
        if superpose_first:
            self.superpose(alignments)
        for i in range(num - 1):
            for j in range(i + 1, num):
                name_1, name_2 = self.structures[i].name, self.structures[j].name
                aln_1 = alignments[name_1]
                aln_2 = alignments[name_2]
                common_coords_1, common_coords_2 = get_common_coordinates(
                    self.structures[i].coordinates,
                    self.structures[j].coordinates,
                    aln_1,
                    aln_2,
                )
                assert common_coords_1.shape[0] > 0
                if not superpose_first:
                    rot, tran = superposition_functions.paired_svd_superpose(
                        common_coords_1, common_coords_2
                    )
                    common_coords_2 = superposition_functions.apply_rotran(
                        common_coords_2, rot, tran
                    )
                pairwise_rmsd_matrix[i, j] = pairwise_rmsd_matrix[
                    j, i
                ] = score_functions.get_rmsd(common_coords_1, common_coords_2)
                pairwise_coverage[i, j] = pairwise_coverage[
                    j, i
                ] = common_coords_1.shape[0] / len(aln_1)
                pairwise_tm[i, j] = pairwise_tm[j, i] = tm_score(
                    common_coords_1,
                    common_coords_2,
                    self.structures[i].length,
                    self.structures[j].length,
                )
        return pairwise_rmsd_matrix, pairwise_coverage, pairwise_tm

    def get_aligned_structures(self, alignment=None):
        if alignment is None:
            alignment = self.alignment
        if type(alignment[self.structures[0].name]) == str:
            alignment = alignment_to_numpy(alignment)
        self.superpose(alignment)
        aligned_structures = np.zeros(
            (len(self.structures), len(alignment[self.structures[0].name]), 3)
        )
        nan_3 = np.zeros(3)
        nan_3[:] = np.nan
        for i in range(len(self.structures)):
            aligned_structures[i] = [
                self.structures[i].coordinates[x] if x != -1 else nan_3
                for x in alignment[self.structures[i].name]
            ]
        return aligned_structures
def trigger_numba_compilation():
    """
    Run this at the beginning of a Caretta run to compile Numba functions
    """
    parameters = {
        "size": 1,
        "gap_open_penalty": 0.0,
        "gap_extend_penalty": 0.0,
        "gamma": 0.03,
    }
    coords_1 = np.zeros((2, 3))
    coords_2 = np.zeros((2, 3))
    tm_score(coords_1, coords_2, 2, 2)
    weights_1 = np.zeros((coords_1.shape[0], 1))
    weights_2 = np.zeros((coords_2.shape[0], 1))
    get_pairwise_alignment(
        coords_1, coords_2, parameters["gamma"], 0, 0, weights_1, weights_2
    )
    superposition_functions.signal_svd_superpose_function(
        coords_1, coords_2, parameters
    )
    distance_matrix = np.random.random((5, 5))
    nj.neighbor_joining(distance_matrix)
    aln_1 = np.array([0, -1, 1])
    aln_2 = np.array([0, 1, -1])
    get_common_coordinates(coords_1, coords_2, aln_1, aln_2)
    get_mean_coords(aln_1, coords_1, aln_2, coords_2)




######## Beginning ########
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help='pdb files list with path included', default='./pdbfile.list')
parser.add_argument('-o', type=str, help='matrix output path', default='./')
args = parser.parse_args()

## collect pdb files
scopepathfile = args.i
pdbs = [line.strip() for line in open(scopepathfile, 'r')]

protein_names = [((pdb.split("/")[-1]).split("."))[0] for pdb in pdbs]

## import pdb files --> change to structure database
msa_class = StructureMultiple.from_pdb_files(scopepathfile)

start_runtime = time()
## run pairwise shape alignment
StructureMultiple.make_pairwise_shape_matrix(msa_class)
matrix = msa_class.pairwise_distance_matrix
end_runtime = time()

print(f"Runtime: {(end_runtime - start_runtime):.2f} seconds")

## convert matrix to 3 column list (query, target, value)
outputfile = args.o
output = open(outputfile, 'w')

for i in range(len(protein_names)):
    query = protein_names[i]
    for j in range(len(protein_names)):
        target = protein_names[j]
        value = str(matrix[i][j])
        output.write(query + "\t" + target + "\t" + value + "\n")
output.close()

