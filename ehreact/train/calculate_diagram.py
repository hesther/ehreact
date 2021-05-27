"""
calculate_diagram.py
Training Hasse diagrams.
"""

import time
import numpy as np
import pickle

from .hasse import extended_hasse
from ehreact.helpers.rdkit import (
    read_in_smiles_unique,
    read_in_reactions_unique,
    preprocess_seeds,
    make_fragment_indices,
    morgan_bit_fp_from_mol,
    tanimoto_from_fp,
    find_matching_atoms,
    canonicalize,
    make_mol,
    make_mol_no_sanit_h,
)
from ehreact.diagram.diagram import sanity_check
from ehreact.diagram.plot_hasse import plot_hasse
from ehreact.helpers.utils import findsubsets
from ehreact.helpers.transition_state import compute_aam_with_h


def calculate_diagram(
    smiles,
    verbose=False,
    quiet=True,
    compute_aam=False,
    save_path=None,
    save_plot=None,
    train_mode="transition_state",
    seed=[],
    no_props=False,
    plot_only_branches=False,
    temp_dir_img=None,
):
    """
    Computes a Hasse diagram of a list of reaction or molecule smiles.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES or reaction SMILES.
    verbose: bool, default False
        Whether to print additional information.
    quiet: bool, default True
        Whether to silence all output.
    compute_aam: bool, default False
        Whether to compute atom-mappings for reactions.
    save_path: str, default  None
        File to which diagram is saved.
    save_plot:  str, default  None
        File to which save image of diagram.
    train_mode: Literal["single_reactant","transition_state"], default "transition_state"
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.
    seed: List[str], default []
        List of SMILES seeds for the reactant algorithm, usually a single seed is given.
    no_props: bool, default False
        Do not compute any properties, just output the diagram.
    plot_only_branches: bool, default False
        Plot only substructures that branch off.
    temp_dir_img: str, default None
        Directory to save temporary image files

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The Hasse diagram of the input list of molecules/reactions.
    """

    # Create Hasse diagram
    if not quiet:
        print("=" * 40)
        print("Calculating diagram")
    start1 = time.time()
    if train_mode == "single_reactant":
        d, smiles_dict = calculate_diagram_single_reactant(smiles, seed, verbose, quiet)
    elif train_mode == "transition_state":
        d, smiles_dict = calculate_diagram_transition_state(
            smiles, verbose, quiet, compute_aam
        )
    else:
        raise ValueError("Unknown option used for train_mode.")
    d.mode = train_mode
    sanity_check(d)
    end1 = time.time()

    # Fill additional information into diagram
    if not quiet:
        print("=" * 40)
        print("Create additional diagram information")
    start2 = time.time()
    # Need to calculate the lowest template for each node for correct plotting of diagram
    for n in d.nodes:
        d.nodes[n].lowest_template = find_lowest_template(d.nodes[n], d)
    # Calculate all other properties only if not no_props:
    if not no_props:
        fill_information(d, train_mode, verbose, smiles_dict)
    end2 = time.time()

    # Save diagram
    start3 = time.time()
    if save_path is not None:
        if not quiet:
            print("=" * 40)
            print("Saving model")
        with open(save_path, "wb") as handle:
            pickle.dump(d, handle)
        if not quiet:
            print("Saved model to", save_path)
    end3 = time.time()

    # Plot diagram
    start4 = time.time()
    if save_plot is not None:
        if not quiet:
            print("=" * 40)
            print("Plotting model")
        if temp_dir_img is None:
            raise ValueError(
                "Must supply a directory for saving temporary images files."
            )
        plot_hasse(d, temp_dir_img, save_plot, plot_only_branches)
        if not quiet:
            print("Plotted diagram in png format to " + save_plot)
    end4 = time.time()

    # Print timing
    if not quiet:
        print("=" * 40)
        print("Finished")
        print("Calculating diagram -- time:{}s".format(end1 - start1))
        print("Create additional diagram information -- time:{}s".format(end2 - start2))
        print("Saving model -- time:{}s".format(end3 - start3))
        print("Plotting model -- time:{}s".format(end4 - start4))

    return d


def calculate_diagram_single_reactant(smiles, seed_list, verbose, quiet):
    """
    Computes a Hasse diagram of a list of molecule smiles.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES.
    seed_list: List[str]
        List of SMILES seeds.
    verbose: bool
        Whether to print additional information.
    quiet: bool
        Whether to silence all output.

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The Hasse diagram of the input list of molecules.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    """

    # Create canonical smiles and mol objects:
    smiles, smiles_dict, skipped = read_in_smiles_unique(smiles)
    if not quiet and len(skipped) != 0:
        print("Skipped duplicates:", skipped)
    if verbose:
        print(
            "Input:",
            len(smiles) + len(skipped),
            "(of which",
            len(skipped),
            "were skipped due to duplicates).",
        )
        print("Canonicalized smiles used for training:")
        for smi in smiles:
            print(
                "Without stereochemistry:",
                smi,
                " with stereochemistry:",
                smiles_dict[smi]["canonical_smi"],
            )

    # Create seeds:
    if verbose and len(seed_list) == 0:
        print(
            "No seed given, find maximum common substructure of all molecules instead."
        )
    seeds, rule_dict, num_smiles_seed = preprocess_seeds(seed_list, smiles, smiles_dict)

    if verbose:
        print("Seeds", seeds)

    if np.sum(num_smiles_seed) < len(smiles):
        raise ValueError("Not all molecules fit a seed. Check given seeds.")
    elif np.sum(num_smiles_seed) > len(smiles):
        raise ValueError("Seeds are not mutually exclusive. Check given seeds.")

    tags_core = {}
    d = extended_hasse(smiles_dict, seeds, rule_dict, tags_core, verbose, quiet)

    return d, smiles_dict


def calculate_diagram_transition_state(smiles, verbose, quiet, compute_aam):
    """
    Computes a Hasse diagram of a list of reaction smiles.

    Parameters
    ----------
    smiles: List[str]
        List of reaction SMILES.
    verbose: bool, default False
        Whether to print additional information.
    quiet: bool
        Whether to silence all output.
    compute_aam: bool
        Whether to compute atom-mappings for reactions.

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The Hasse diagram of the input list of reactions.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    """

    # Create transition states
    if compute_aam:
        if not quiet:
            print("Calculating atom-mappings")
        mapped_smiles = []
        for smi in smiles:
            mapped_smiles.append(compute_aam_with_h(smi))
            if verbose:
                print(smi, mapped_smiles[-1])
        smiles = mapped_smiles
        if not quiet:
            print("Finished calculating atom-mappings")
    (
        seeds,
        rule_dict,
        smiles,
        smiles_dict,
        skipped,
        tags_core,
    ) = read_in_reactions_unique(smiles)
    if not quiet and len(skipped) != 0:
        print("Skipped duplicates:", skipped)
    if verbose:
        print(
            "Input:",
            len(smiles) + len(skipped),
            "(of which",
            len(skipped),
            "were skipped due to duplicates).",
        )
        print("Canonicalized smiles used for training:")
        for smi in smiles:
            print(
                "Without stereochemistry:",
                smi,
                " with stereochemistry:",
                smiles_dict[smi]["canonical_smi"],
            )

    d = extended_hasse(smiles_dict, seeds, rule_dict, tags_core, verbose, quiet)

    return d, smiles_dict


def find_lowest_template(curr_node, d):
    """Function to find the lowest (most general) substructure/reaction rule in the tree.

    Parameters
    ----------
    curr_node: ehreact.diagram.diagram.Node
         Node for which to find the lowest template.
    d: ehreact.diagram.diagram.Diagram
         Hasse diagram.

    Returns
    -------
    lowest_template: str
        Name of the lowest template.
    """

    if curr_node.key == "":
        return ""
    while curr_node.edges_to_parent[0].parent_node.key != "":
        curr_node = curr_node.edges_to_parent[0].parent_node
    lowest_template = curr_node.key
    return lowest_template


def fill_information(d, train_mode, verbose, smiles_dict):
    """Function to fill topology information and fingerprints into a Hasse diagram (alters diagram in-place).

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
         Hasse diagram.
    train_mode: Literal["single_reactant","transition_state"]
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.
    verbose: bool, default False
        Whether to print additional information.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    """

    # Make a leaf and non_leaf list, since they need different information to be created
    leafs = []
    non_leafs = []
    for n in d.nodes:
        if d.nodes[n].is_leaf:
            leafs.append(d.nodes[n].key)
        elif d.nodes[n].key != "":
            non_leafs.append(d.nodes[n].key)

    # For all leafs, calculate the distance to every non-leaf
    if verbose:
        print("... calculating minimum distance to leaf nodes")
    for leaf in leafs:
        d.nodes[leaf].min_dist_leaf = 0
        child = d.nodes[leaf]
        while child.key != "":
            parent = child.edges_to_parent[0].parent_node
            parent.min_dist_leaf = min(parent.min_dist_leaf, child.min_dist_leaf + 1)
            parent.all_leafs.append(leaf)
            child = parent
    if verbose:
        for non_leaf in non_leafs:
            print("...... distance", d.nodes[non_leaf].min_dist_leaf, "for", non_leaf)

    # Calculate Morgan fingerprints for each leaf
    if verbose:
        print("Calculating fingerprints")
    for leaf in leafs:
        mol = smiles_dict[leaf]["mol"]
        mol_no_stereo = smiles_dict[leaf]["mol_no_stereo"]
        if train_mode == "single_reactant":
            fp = morgan_bit_fp_from_mol(mol_no_stereo)
            fp_stereo = [morgan_bit_fp_from_mol(item, chirality=True) for item in mol]
            d.nodes[leaf].fp_reac = fp
            d.nodes[leaf].fp_reac_stereo = fp_stereo
        elif train_mode == "transition_state":
            fp_reac = morgan_bit_fp_from_mol(mol_no_stereo[0])
            fp_reac_stereo = [
                morgan_bit_fp_from_mol(item[0], chirality=True) for item in mol
            ]
            d.nodes[leaf].fp_reac = fp_reac
            d.nodes[leaf].fp_reac_stereo = fp_reac_stereo
            fp_prod = morgan_bit_fp_from_mol(mol_no_stereo[1])
            fp_prod_stereo = [
                morgan_bit_fp_from_mol(item[1], chirality=True) for item in mol
            ]
            d.nodes[leaf].fp_prod = fp_prod
            d.nodes[leaf].fp_prod_stereo = fp_prod_stereo

    # Calculate overall diversity of diagram:
    if verbose:
        print("Calculating overall diversity")
    div_mean_reac, div_mean_prod = calculate_diversity(d, "", False, train_mode)
    d.nodes[""].diversity_reac = div_mean_reac
    d.nodes[""].diversity_prod = div_mean_prod
    div_mean_reac, div_mean_prod = calculate_diversity(d, "", True, train_mode)
    d.nodes[""].diversity_reac_stereo = div_mean_reac
    d.nodes[""].diversity_prod_stereo = div_mean_prod

    # Write fragment list to root
    if train_mode == "transition_state":
        if verbose:
            print("Creating fragment list")
        write_fragment_list_to_root(d, train_mode, verbose, smiles_dict)


def write_fragment_list_to_root(d, train_mode, verbose, smiles_dict):
    """Function to calculate a list of reactant rule fragments (only atoms in reaction center). This in needed
    to transform inputted molecules to their corresponding transition state. Save to the root node (in-place).

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
         Hasse diagram.
    train_mode: Literal["single_reactant","transition_state"]
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.
    verbose: bool, default False
        Whether to print additional information.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    """

    if verbose:
        print("Calculating reactant fragments for each minimal template")
    fragment_dict_reac = {}
    minimal_templates_edges = d.nodes[""].edges_to_child
    for edge in minimal_templates_edges:
        minimal_template = edge.child_node
        if verbose:
            print("template", minimal_template)
        fragment_dict_reac[str(minimal_template)] = {}
        rule_reac = d.nodes[str(minimal_template)].rule_reac
        rule_reac_fragment_indices, rule_reac_fragment_mols = make_fragment_indices(
            rule_reac
        )
        for i in range(len(rule_reac_fragment_indices)):
            fragment_dict_reac[str(minimal_template)][rule_reac_fragment_indices[i]] = {
                "rule_fragment": rule_reac_fragment_mols[i],
                "molecules": [],
            }
        for leaf in d.nodes[str(minimal_template)].all_leafs:
            if verbose:
                print("... leaf", leaf)
            tags_core = d.nodes[leaf].tags_core
            smi_reac = leaf.split(">")[0]
            mol_reac = make_mol_no_sanit_h(smi_reac)
            matches_reac = find_matching_atoms(
                train_mode, mol_reac, rule_reac, None, tags_core
            )
            if len(matches_reac) != 1:
                print(
                    "WARNING: Expected exactly one match, but got "
                    + str(len(matches_reac))
                    + " instead."
                )
                print("tags:", tags_core)
                print("leaf:", leaf)
                print("minimal_template:", minimal_template)
            if len(matches_reac) == 0:
                raise ValueError("Zero matches, something went wrong.")
            matches_reac = matches_reac[0]
            if len(smi_reac.split(".")) != len(rule_reac_fragment_indices):
                print(
                    "Intramolecular reaction detected. This is currently not tested to work correctly"
                )
                print("this might lead to erraneous results.")
                print(leaf)
                print(smi_reac)
                print(rule_reac_fragment_indices)
                # TODO implement intramolecular version, write exception if rule_atoms_in_fragment not found
                # (since it is a superset of the combination of available fragments then,
                # need to iterate through those!
            for smi_fragment in smi_reac.split("."):
                fragment = make_mol(smi_fragment)
                fragment_matches = mol_reac.GetSubstructMatches(fragment)
                fragment_match = []
                for match in fragment_matches:
                    correct_match = True
                    for i, x in enumerate(match):
                        if fragment.GetAtomWithIdx(i).HasProp("molAtomMapNumber"):
                            if fragment.GetAtomWithIdx(i).GetProp(
                                "molAtomMapNumber"
                            ) != mol_reac.GetAtoms()[x].GetProp("molAtomMapNumber"):
                                correct_match = False
                    if correct_match:
                        fragment_match.append(match)
                if len(fragment_match) != 1:
                    raise ValueError(
                        "Expected exactly one match, but got",
                        len(fragment_match),
                        "instead. Exit program",
                    )
                fragment_match = fragment_match[0]
                rule_atoms_in_fragment = tuple(
                    [i for i, x in enumerate(matches_reac) if x in fragment_match]
                )
                for atom in fragment.GetAtoms():
                    atom.SetAtomMapNum(0)
                canonical_smi_fragment = canonicalize(fragment)
                if (
                    rule_atoms_in_fragment
                    not in fragment_dict_reac[str(minimal_template)].keys()
                ):
                    # Could be in intramolecular reaction, correct for it here:
                    # TODO, add to molecules for each of the subsets that apply.
                    # else this is an error:
                    raise ValueError(
                        "An error occured: Found rule fragment doesn't match available fragments"
                    )
                if (
                    canonical_smi_fragment
                    not in fragment_dict_reac[str(minimal_template)][
                        rule_atoms_in_fragment
                    ]["molecules"]
                ):
                    fragment_dict_reac[str(minimal_template)][rule_atoms_in_fragment][
                        "molecules"
                    ].append(canonical_smi_fragment)

    d.nodes[""].fragment_reac = fragment_dict_reac


def calculate_diversity(d, node, stereo, train_mode):
    """Calculates diversity within a branch or tree.

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
         Hasse diagram.
    node: ehreact.diagram.diagram.Node
         Node for which to calculate diversity, all leaf nodes attached to this node by an arbitrary number of
         edges toward children are taken into account.
    stereo: bool
         Whether to include stereochemistry in fingerprints
    train_mode: Literal["single_reactant", "transition_state"]
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.

    Returns
    -------
    mean_div_reac: float
        Mean pair similarity of reactants.
    mean_div_prod: float
        Mean pair similarity of products.
    """

    diversity_reac = []
    diversity_prod = []
    if len(d.nodes[node].all_leafs) > 1:
        for combi in findsubsets(d.nodes[node].all_leafs, 2):
            if not stereo:
                fp1 = d.nodes[combi[0]].fp_reac
                fp2 = d.nodes[combi[1]].fp_reac
                diversity_reac.append(tanimoto_from_fp(fp1, fp2))
                if train_mode == "transition_state":
                    fp1 = d.nodes[combi[0]].fp_prod
                    fp2 = d.nodes[combi[1]].fp_prod
                    diversity_prod.append(tanimoto_from_fp(fp1, fp2))
                else:
                    diversity_prod.append(1)
            else:
                best_sim = 0
                for fp1 in d.nodes[combi[0]].fp_reac_stereo:
                    for fp2 in d.nodes[combi[1]].fp_reac_stereo:
                        sim = tanimoto_from_fp(fp1, fp2)
                        if sim > best_sim:
                            best_sim = sim
                diversity_reac.append(best_sim)
                if train_mode == "transition_state":
                    best_sim = 0
                    for fp1 in d.nodes[combi[0]].fp_prod_stereo:
                        for fp2 in d.nodes[combi[1]].fp_prod_stereo:
                            sim = tanimoto_from_fp(fp1, fp2)
                            if sim > best_sim:
                                best_sim = sim
                    diversity_prod.append(best_sim)
                else:
                    diversity_prod.append(1)
    else:
        diversity_reac.append(1)
        diversity_prod.append(1)

    mean_div_reac = np.mean(diversity_reac)
    mean_div_prod = np.mean(diversity_prod)

    return mean_div_reac, mean_div_prod
