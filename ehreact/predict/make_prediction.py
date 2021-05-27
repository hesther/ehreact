"""
Score on Hasse diagram.
"""

import time
from ehreact.helpers.rdkit import (
    read_in_smiles,
    read_in_reactions,
    mol_to_rxn_smiles,
    morgan_bit_fp_from_mol,
    tanimoto_from_fp,
)
from .score import default_score_generator
import numpy as np


def make_prediction(
    smiles,
    d,
    params=None,
    verbose=False,
    quiet=True,
    compute_aam=False,
    predict_mode="transition_state",
    stereochemistry=False,
):
    """
    Computes a Hasse diagram of a list of reaction or molecule smiles.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES or reaction SMILES to score.
    d: ehreact.diagram.diagram.Diagram
        Diagram to calculate scores on.
    params: dict, default None
        Dictionary of scoring hyperparameters.
    verbose: bool, default False
        Whether to print additional information.
    quiet: bool, default True
        Whether to silence all output.
    compute_aam: bool, default False
        Whether to compute atom-mappings for reactions.
    predict_mode: Literal["single_reactant","multi_reactant","transition_state"], default "transition_state"
        Prediction mode, either transition states extracted from reaction smiles or single/multi reactants from smiles.
        If single reactants, reaction partners are added automatically and all possible reaction products enumerated,
        for multi_reactants no reaction partners are added but different reaction outcomes are taken in to account.
    stereochemistry: bool, default False
        Whether to use stereochemistry for scoring.

    Returns
    -------
    scores: List[float]
        List of scores.
    combination: List[str]
        List of multiple reactants.
    current_smiles: List[str]
        List of SMILES strings.
    belongs_to: List[int]
        List of integer index characterizing to which initial smiles the reaction belongs to.
    raw_scores: List[dict]
        List of dictionaries of raw scores.
    """

    # Initialize dictionaries and lists
    if params is None:
        params={"use_prod": True,"coeff_ss": 1.0,"coeff_sp": -1.0,"coeff_sm": 1.0,"coeff_sl": -0.1,"cap_sp": 0.8,"cap_sl": 5.0}

    if not quiet:
        print("Parameters for scoring:")
        print(params)

    start1 = time.time()

    if not quiet:
        print("=" * 40)
        print("Preprocessing data")

    tags_core = {}
    combination = [""] * len(smiles)
    belongs_to = list(range(len(smiles)))

    # Preprocess for single_reactant training mode
    if d.mode == "single_reactant":
        if predict_mode != "single_reactant":
            raise ValueError(
                "In train_mode='single_reactant' it is required to also use predict_mode='single_reactant'."
            )
        current_smiles, smiles_dict = read_in_smiles(smiles)
        current_mols = [smiles_dict[smi]["mol_diagram"] for smi in current_smiles]
        for i in range(len(current_smiles)):
            tags_core[current_smiles[i]] = []

    # Preprocess for transition_state training mode
    elif d.mode == "transition_state":
        if predict_mode == "single_reactant" or predict_mode == "multi_reactant":
            if stereochemistry:
                raise ValueError(
                    """Use of stereochemistry in scores is currently only available for
                    single_reactant/single_reactant or transition_state/transition_state mode."""
                )
            initial_smiles_no_stereo, initial_smiles_dict = read_in_smiles(smiles)
            initial_mols = [
                initial_smiles_dict[smi]["mol_diagram"]
                for smi in initial_smiles_no_stereo
            ]

            # TODO intramolecular reactions
            current_smiles, belongs_to, combination = mol_to_rxn_smiles(
                initial_smiles_no_stereo, initial_mols, d, predict_mode, verbose
            )

            if verbose:
                print(
                    "Inputting smiles string to transition state routine",
                    current_smiles,
                )
            current_smiles, smiles_dict, tags_core = read_in_reactions(current_smiles)
            current_mols = [smiles_dict[smi]["mol_diagram"] for smi in current_smiles]

        elif predict_mode == "transition_state":

            current_smiles, smiles_dict, tags_core = read_in_reactions(smiles)
            current_mols = [smiles_dict[smi]["mol_diagram"] for smi in current_smiles]

    end1 = time.time()

    # Prediction on diagram
    if not quiet:
        print("=" * 40)
        print("Prediction on diagram")

    start2 = time.time()

    scores = []
    raw_scores = []
    for i in range(len(current_smiles)):
        if verbose:
            print("Searching template for", current_smiles[i])
        highest_template = find_highest_template(
            curr_node=d.nodes[""], mol=current_mols[i], d=d, highest_template=None
        )
        if highest_template != "":
            if verbose:
                print(
                    "Highest template:",
                    d.nodes[highest_template].key,
                    "with",
                    d.nodes[highest_template].min_dist_leaf,
                    "edges to nearest leaf and",
                    len(d.nodes[highest_template].all_leafs),
                    "leaf nodes",
                )
            score, raw_score = calculate_score(
                current_smiles[i],
                current_mols[i],
                highest_template,
                d,
                verbose,
                predict_mode,
                params,
                smiles_dict,
                stereochemistry,
                tags_core=tags_core[current_smiles[i]],
            )
        else:
            if verbose:
                print("No template found, setting score to 0.")
            score, raw_score = 0, {}
        scores.append(score)
        raw_scores.append(raw_score)

    if not quiet:
        print("=" * 40)
        print("Scores:")
        for i in range(len(current_smiles)):
            print(
                smiles[belongs_to[i]],
                combination[i],
                current_smiles[i],
                ": Score:",
                scores[i],
            )

    end2 = time.time()

    # Print timing
    if not quiet:
        print("=" * 40)
        print("Finished")
        print("Preprocessing -- time:{}s".format(end1 - start1))
        print("Prediction -- time:{}s".format(end2 - start2))

    return scores, combination, current_smiles, belongs_to, raw_scores


def find_highest_template(curr_node, mol, d, highest_template):
    """
    Recursive function to find the highest matching substructure/reaction rule of a molecule/reaction.

    Parameters
    ----------
    curr_node: str
        SMILES string of current node.
    mol: rdkit.Chem.Mol
        RDKit molecule.
    d: ehreact.diagram.diagram.Diagram
        Hasse diagram.
    highest_template: str
        Current highest matching substructure/reaction rule.

    Returns
    -------
    highest_template: str
        highest matching (most specific) template in tree
    """

    if curr_node.key != "":
        rule = curr_node.rule
        match = mol.GetSubstructMatch(rule)
    else:
        match = True
    if match or curr_node.key == "":
        highest_template = curr_node.key
        for e in curr_node.edges_to_child:
            n = e.child_node
            if not n.is_leaf:
                highest_template = find_highest_template(n, mol, d, highest_template)
    return highest_template


def calculate_score(
    smi,
    mol,
    highest_template,
    d,
    verbose,
    predict_mode,
    params,
    smiles_dict,
    stereochemistry,
    tags_core=None,
):
    """
    Function to score a query molecule given the diagram, highest template, QM information.

    Parameters
    ----------
    smi: str
        SMILES string of the query molecule.
    mol: rdkit.Chem.Mol
        RDKit molecule.
    highest_template: str
        key of the highest matching template.
    d: ehreact.diagram.diagram.Diagram
        Hasse diagram.
    verbose: bool, default False
        Whether to print additional information.
    predict_mode: Literal["single_reactant","multi_reactant","transition_state"], default "transition_state"
        Prediction mode, either transition states extracted from reaction smiles or single/multi reactants from smiles.
        If single reactants, reaction partners are added automatically and all possible reaction products enumerated,
        for multi_reactants no reaction partners are added but different reaction outcomes are taken in to account.
     params: dict
        Dictionary of scoring hyperparameters.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    stereochemistry: bool
        Whether to use stereochemistry for scoring.
    tags_core: dict, default None
        Atom mapping numbers of atom in the reaction core. Used to verify a correct match.

    Returns
    -------
    best_score: float
        Score.
    rawnumber: dict
        Dictionary of raw scores.
    """

    if verbose:
        print("Scoring", smi, "at template", highest_template)

    overall_div_mean_reac = d.nodes[""].diversity_reac
    overall_div_mean_prod = d.nodes[""].diversity_prod

    chemical_scores = []
    chemical_scores_prod = []
    chemical_scores_overall = []
    chemical_scores_overall_prod = []

    if predict_mode == "single_reactant" and d.mode == "single_reactant":
        if stereochemistry:
            fp_reac = morgan_bit_fp_from_mol(smiles_dict[smi]["mol"][0], chirality=True)
        else:
            fp_reac = morgan_bit_fp_from_mol(smiles_dict[smi]["mol_no_stereo"])

        for leaf in d.nodes[highest_template].all_leafs:
            if stereochemistry:
                fp_reac_ref = d.nodes[leaf].fp_reac_stereo
                chemical_score = np.max(
                    [tanimoto_from_fp(fp_reac, ref) for ref in fp_reac_ref]
                )
            else:
                fp_reac_ref = d.nodes[leaf].fp_reac
                chemical_score = tanimoto_from_fp(fp_reac, fp_reac_ref)
            chemical_scores.append(chemical_score)

        for leaf in d.nodes[""].all_leafs:
            if stereochemistry:
                fp_reac_ref = d.nodes[leaf].fp_reac_stereo
                chemical_score = np.max(
                    [tanimoto_from_fp(fp_reac, ref) for ref in fp_reac_ref]
                )
            else:
                fp_reac_ref = d.nodes[leaf].fp_reac
                chemical_score = tanimoto_from_fp(fp_reac, fp_reac_ref)
            chemical_scores_overall.append(chemical_score)

    elif d.mode == "transition_state":
        if stereochemistry:
            fp_reac = morgan_bit_fp_from_mol(
                smiles_dict[smi]["mol"][0][0], chirality=True
            )
        else:
            fp_reac = morgan_bit_fp_from_mol(smiles_dict[smi]["mol_no_stereo"][0])
        if params["use_prod"]:
            if stereochemistry:
                fp_prod = morgan_bit_fp_from_mol(
                    smiles_dict[smi]["mol"][0][1], chirality=True
                )
            else:
                fp_prod = morgan_bit_fp_from_mol(smiles_dict[smi]["mol_no_stereo"][1])

        for leaf in d.nodes[highest_template].all_leafs:
            if stereochemistry:
                fp_reac_ref = d.nodes[leaf].fp_reac_stereo
                fp_prod_ref = d.nodes[leaf].fp_prod_stereo
            else:
                fp_reac_ref = d.nodes[leaf].fp_reac
                fp_prod_ref = d.nodes[leaf].fp_prod
            if stereochemistry:
                chemical_score = np.max(
                    [tanimoto_from_fp(fp_reac, ref) for ref in fp_reac_ref]
                )
            else:
                chemical_score = tanimoto_from_fp(fp_reac, fp_reac_ref)
            chemical_scores.append(chemical_score)
            if params["use_prod"]:
                if stereochemistry:
                    chemical_score_prod = np.max(
                        [tanimoto_from_fp(fp_prod, ref) for ref in fp_prod_ref]
                    )
                else:
                    chemical_score_prod = tanimoto_from_fp(fp_prod, fp_prod_ref)
                chemical_scores_prod.append(chemical_score_prod)

        for leaf in d.nodes[""].all_leafs:
            if stereochemistry:
                fp_reac_ref = d.nodes[leaf].fp_reac_stereo
                chemical_score = np.max(
                    [tanimoto_from_fp(fp_reac, ref) for ref in fp_reac_ref]
                )
            else:
                fp_reac_ref = d.nodes[leaf].fp_reac
                chemical_score = tanimoto_from_fp(fp_reac, fp_reac_ref)
            chemical_scores_overall.append(chemical_score)
            if params["use_prod"]:
                if stereochemistry:
                    fp_prod_ref = d.nodes[leaf].fp_prod_stereo
                    chemical_score_prod = np.max(
                        [tanimoto_from_fp(fp_prod, ref) for ref in fp_prod_ref]
                    )
                else:
                    fp_prod_ref = d.nodes[leaf].fp_prod
                    chemical_score_prod = tanimoto_from_fp(fp_prod, fp_prod_ref)
                chemical_scores_overall_prod.append(chemical_score_prod)

    min_dist_leaf = d.nodes[highest_template].min_dist_leaf

    rawnumber = {
        "chemical_scores": chemical_scores,
        "chemical_scores_prod": chemical_scores_prod,
        "chemical_scores_overall": chemical_scores_overall,
        "chemical_scores_overall_prod": chemical_scores_overall_prod,
        "overall_div_mean_reac": overall_div_mean_reac,
        "overall_div_mean_prod": overall_div_mean_prod,
        "min_dist_leaf": min_dist_leaf,
        "train_mode": d.mode,
    }

    best_score = default_score_generator(rawnumber, params, verbose)

    return best_score, rawnumber
