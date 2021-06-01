from ehreact.helpers.rdkit import moltosmiles_transition, mol_transition
from ehreact.diagram.diagram import Diagram
from copy import deepcopy
from rdkit import Chem
from ehreact.helpers.utils import findsubsets
import numpy as np

def extended_hasse(smiles_dict, seeds, rule_dict, tags_core, verbose, quiet):
    """
    Create an extended Hasse diagram.

    Parameters
    ----------
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    seeds: List[str]
        List of SMILES seeds.
    rule_dict: dict
        A dictionary of all minimal templates of all seeds.
    tags_core: dict
        A dictionary of the atom map numbers for each minimal template. Empty dictionary for single reactant mode.
    verbose: bool
        Whether to print additional information.
    quiet: bool
        Whether to silence all output.

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The Hasse diagram of the input list of molecules.
    """

    d = Diagram()  # create empty diagram
    for seed in seeds:
        change_dict = rule_dict[seed]["change_dict"]
        smiles = rule_dict[seed]["smiles"]
        mols = [smiles_dict[smi]["mol_diagram"] for smi in smiles]
        if len(smiles) != 0:
            rule = rule_dict[seed]["rule"]
            rule_smiles = moltosmiles_transition(rule, change_dict)
            if not quiet:
                print("Calculating diagram for seed", seed)
            # Add root node:
            d.add_node(
                key_node=rule_smiles,
                rule=rule,
                is_leaf=False,
                key_parent_node="",
                quiet=quiet,
            )
            d.nodes[rule_smiles].rule_reac = mol_transition(rule, change_dict["reac"])
            d.nodes[rule_smiles].rule_prod = mol_transition(rule, change_dict["prod"])
            d.nodes[""].change_dict[rule_smiles] = change_dict

            # Start iteration through diagram:
            iterate_algorithm(
                d,
                mols,
                smiles,
                rule,
                rule_smiles,
                verbose,
                quiet,
                change_dict,
                tags_core,
            )

    return d


def iterate_algorithm(
    d, mols, smiles, rule, rule_smiles, verbose, quiet, change_dict, tags_core
):
    """
    Iterate atom extension algorithm by looking for possible atoms to add to, find best combination,
    create a new template and attach it to the current Hasse diagram (edits to diagram are in-place).

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
        The current Hasse diagram.
    mols: List[rdkit.Chem.Mol]
        List of RDKit molecule objects of the molecules or pseudomolecules which are considered in the current branch.
    smiles: List[str]
        List of SMILES strings of the molecules or pseudomolecules which are considered in the current branch.
    rule: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    rule_smiles: str
        Name of the current template.
    verbose: bool
        Whether to print additional information.
    quiet: bool
        Whether to silence all output.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    tags_core: dict
        A dictionary of the atom map numbers for each minimal template. Empty dictionary for single reactant mode.
    """

    if verbose:
        print("...... NEW ITERATION with mols", smiles)
    new_rules = enlarge_rule(rule, mols, smiles, change_dict, verbose)
    for new_rule_smiles in new_rules.keys():
        new_rule = new_rules[new_rule_smiles]["rule"]
        new_mols = new_rules[new_rule_smiles]["mols"]
        new_smiles = new_rules[new_rule_smiles]["smiles"]
        if len(new_mols) == 1 and new_rule.GetNumAtoms() == new_mols[0].GetNumAtoms():
            # leaf node detected, no further iteration in this branch
            d.add_node(
                key_node=new_smiles[0],
                rule=new_rule,
                is_leaf=True,
                key_parent_node=rule_smiles,
                quiet=quiet,
            )
            if tags_core != {}:  # transition_state_mode:
                d.nodes[new_smiles[0]].tags_core = tags_core[new_smiles[0]]
        else:
            d.add_node(
                key_node=new_rule_smiles,
                rule=new_rule,
                is_leaf=False,
                key_parent_node=rule_smiles,
                quiet=quiet,
            )
            if tags_core != {}:  # transition_state_mode:
                d.nodes[new_rule_smiles].rule_reac = mol_transition(
                    new_rule, change_dict["reac"]
                )
                d.nodes[new_rule_smiles].rule_prod = mol_transition(
                    new_rule, change_dict["prod"]
                )
            iterate_algorithm(
                d,
                new_mols,
                new_smiles,
                new_rule,
                new_rule_smiles,
                verbose,
                quiet,
                change_dict,
                tags_core,
            )


def enlarge_rule(rule, mols, smiles, change_dict, verbose):
    """
    Function to look for possible atoms to add to current template, find best combination,
    and create a new template.

    Parameters
    ----------
    rule: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    mols: List[rdkit.Chem.Mol]
        List of RDKit molecule objects of the molecules or pseudomolecules which are considered in the current branch.
    smiles: List[str]
        List of SMILES strings of the molecules or pseudomolecules which are considered in the current branch.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    verbose: bool
        Whether to print additional information.

    Returns
    -------
    new_rules: dict
        A dictionary of new templates (childs to the current node), might be one or multiple templates.
    """

    new_rules = {}
    possible_extensions = get_possible_extensions(rule, mols, smiles)
    max_possible_rules = get_max_possible_rules(
        possible_extensions, rule, mols, change_dict
    )
    if verbose:
        print("...... Possible extensions: ", possible_extensions)
        print("...... Possible max. rules: ", max_possible_rules)

    if len(mols) == 1:  # Simply extend all atoms
        match = list(possible_extensions[smiles[0]].keys())[0]
        if verbose:
            print("...... Single molecule in branch, extend at all atoms")
        new_rule = list(max_possible_rules[smiles[0]][match].values())[0]
        new_rule_smiles = moltosmiles_transition(new_rule, change_dict)
        new_rules.update(
            {new_rule_smiles: {"rule": new_rule, "mols": mols, "smiles": smiles}}
        )

    else:
        # use molecule with least number of possible extensions to speed up matching:
        min_num_ext = 999
        min_num_matches = 999
        min_match = None
        for i in range(len(mols)):
            matches = list(possible_extensions[smiles[i]].keys())
            num_matches = len(matches)
            for match in matches:
                num_ext = len(possible_extensions[smiles[i]][match].keys())
                if num_ext < min_num_ext:
                    min_num_idx = i
                    min_num_ext = num_ext
                    min_num_matches = num_matches
                    min_match = match
                elif num_ext == min_num_ext:
                    if num_matches < min_num_matches:
                        min_num_idx = i
                        min_num_ext = num_ext
                        min_num_matches = num_matches
                        min_match = match
        if verbose:
            print(
                "...... Taking molecule",
                smiles[min_num_idx],
                "as pivot, with",
                min_num_matches,
                "matches.",
            )

        if num_matches == 1:
            new_rules = extend_by_single_match(
                possible_extensions,
                min_match,
                min_num_idx,
                mols,
                smiles,
                rule,
                change_dict,
                verbose,
                max_possible_rules,
            )
        else:
            if verbose:
                print(
                    "...... More than one match in pivot => taking the match producing less split branches"
                )
            matches = list(possible_extensions[smiles[min_num_idx]].keys())
            new_rules = extend_by_single_match(
                possible_extensions,
                matches[0],
                min_num_idx,
                mols,
                smiles,
                rule,
                change_dict,
                verbose,
                max_possible_rules,
            )
            if verbose:
                print(new_rules)
            for idx_match in range(1, num_matches):
                try_new_rules = extend_by_single_match(
                    possible_extensions,
                    matches[idx_match],
                    min_num_idx,
                    mols,
                    smiles,
                    rule,
                    change_dict,
                    verbose,
                    max_possible_rules,
                )
                if verbose:
                    print("Try", idx_match + 1, ":", try_new_rules)
                if try_new_rules.keys() != new_rules.keys():
                    if len(try_new_rules.keys()) < len(new_rules.keys()):
                        new_rules = try_new_rules
                        if verbose:
                            print(
                                "...... Updating to new possibility of match", idx_match
                            )
    return new_rules


def check_one_extension(pivot_rule, max_possible_rules):
    """
    Function to iterate over all current molecules/pseudomolecules and check whether each of them
    allows for an extension resulting in the current pivot_rule.

    Parameters
    ----------
    pivot_rule: RDKit.Chem.Mol
        A possible new template
    max_possible_rules: dict
        A dictionary of the matching atoms and possible new templates for all molecules/pseudomolecules.

    Returns
    -------
    bool
        Whether or not the pivot rule has a substructure match with any of the possible extensions.
    """

    for smi in max_possible_rules.keys():
        if not has_one_extension(pivot_rule, max_possible_rules[smi]):
            return False
    return True


def has_one_extension(pivot_rule, extensions_all_matches):
    """
    Function to check whether any current template match and corresponding extension results
    in the current pivot_rule.

    Parameters
    ----------
    pivot_rule: RDKit.Chem.Mol
        A possible new template
    extensions_all_matches: dict
        A dictionary of matching atoms and possible new templates.

    Returns
    -------
    bool
        Whether or not the pivot rule has a substructure match with any of the possible extensions.
    """

    for match in extensions_all_matches.keys():
        max_rule = list(extensions_all_matches[match].values())[0]
        if quick_match(pivot_rule,max_rule): #Only compute matching if plausible (based on number of atoms)
            if max_rule.HasSubstructMatch(pivot_rule):
                return True
    return False


def quick_match(mol_small,mol_large):
    """
    Computes whether mol_small could possibly be a subgraph of mol_large based on atom type counts.

    Parameters
    ----------
    mol_small: RDKit.Chem.Mol
        A molecule.
    mol_large: RDKit.Chem.Mol
        A molecule.

    Returns
    -------
    bool
        Whether a subgraph match is possible based on atom type counts.
    """

    if mol_small.GetNumAtoms() > mol_large.GetNumAtoms():
        return False
    
    bincount_small=np.bincount([atom.GetAtomicNum() for atom in mol_small.GetAtoms()])
    bincount_large=np.bincount([atom.GetAtomicNum() for atom in mol_large.GetAtoms()])

    if bincount_small.shape[0] > bincount_large.shape[0]:
        return False
    
    for i in range(bincount_small.shape[0]):
        if bincount_small[i]>bincount_large[i]:
            return False
    return True


def extend_by_single_match(
    possible_extensions,
    match,
    idx,
    mols,
    smiles,
    rule,
    change_dict,
    verbose,
    max_possible_rules,
):
    """
    Enlarges the current template by selecting the best extensions.

    Parameters
    ----------
    possible_extensions: dict
        A dictionary of possible extensions for each input molecule, containing the matching indices, atom indices of
        atoms to extend, and their possible extension as string.
    match: tuple
        Tuple of the matching atom indices for the pivot molecule.
    idx: str
        Index of pivot molecule.
    mols: List[rdkit.Chem.Mol]
        List of RDKit molecule objects of the molecules or pseudomolecules which are considered in the current branch.
    smiles: List[str]
        List of SMILES strings of the molecules or pseudomolecules which are considered in the current branch.
    rule: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    verbose: bool
        Whether to print additional information.
    max_possible_rules: dict
        Dictionary of possibly extensions, containing the atoms indices of the matching atoms, as well as all
        the possible enlarged templates as string and RDKit molecule.

    Returns
    -------
    new_patterns: dict
        Dictionary of new templates, to be attached to current node as children, including list of molecules and smiles
        belonging to each new template.

    Examples
    --------
    >>> smis = ["CCOCO", "CCOOC"]
    >>> mols = [Chem.MolFromSmiles(smi) for smi in smis]
    >>> rule = Chem.MolFromSmiles("CCO", sanitize=False)
    >>> possible_extensions = {'CCOCO': {(0, 1, 2): {2: 'OC'}}, 'CCOOC': {(0, 1, 2): {2: 'OO'}}}
    >>> change_dict={"reac": {"atom": {}, "bond": {}}, "prod": {"atom": {}, "bond": {}}}
    >>> max_possible_rules=get_max_possible_rules(possible_extensions,rule,mols,change_dict)
    >>> extend_by_single_match(possible_extensions,(0,1,2),0,mols,smis,rule,change_dict,False,max_possible_rules)
    {'COCC': {'rule': <rdkit.Chem.rdchem.Mol at 0x7ffc4d7592b0>,
              'mols': [<rdkit.Chem.rdchem.Mol at 0x7ffc4d75eda0>], 'smiles': ['CCOCO']},
     'CCOO': {'rule': <rdkit.Chem.rdchem.Mol at 0x7ffc4d759f30>,
              'mols': [<rdkit.Chem.rdchem.Mol at 0x7ffc4d75e1c0>], 'smiles': ['CCOOC']}}
    """

    new_patterns = {}
    extensions = possible_extensions[smiles[idx]][match]
    min_extensions_pivot = [
        extend_by_atom(rule, mols[idx], match, [key])
        for key in possible_extensions[smiles[idx]][match].keys()
    ]
    independent_extension = [
        check_one_extension(item, max_possible_rules) for item in min_extensions_pivot
    ]

    if True not in independent_extension:
        # If no extension same for all mols: Extend all and create branches
        if verbose:
            print("...... No common patterns, split to branches")
        all_patt_list = []
        all_patt_list_name = []
        for i in range(len(mols)):
            mol = mols[i]
            new_patt_list = []
            new_patt_list_name = []
            for match in possible_extensions[smiles[i]].keys():
                idx_list = possible_extensions[smiles[i]][match].keys()
                new_patt = extend_by_atom(rule, mol, match, idx_list)
                new_patt_name = moltosmiles_transition(new_patt, change_dict)
                if new_patt_name not in new_patt_list_name:
                    new_patt_list.append(new_patt)
                    new_patt_list_name.append(new_patt_name)
            all_patt_list.append(new_patt_list)
            all_patt_list_name.append(new_patt_list_name)
        for i in range(len(mols)):
            mol = mols[i]
            smi = smiles[i]
            if verbose:
                print("...... Splitting branch for molecule", smi)
            if len(all_patt_list[i]) == 1:
                if verbose:
                    print(
                        "...... Only one match, add template", all_patt_list_name[i][0]
                    )
                new_patt = all_patt_list[i][0]
                if not all_patt_list_name[i][0] in new_patterns:
                    new_patterns.update(
                        {
                            all_patt_list_name[i][0]: {
                                "rule": new_patt,
                                "mols": [
                                    mol,
                                ],
                                "smiles": [
                                    smi,
                                ],
                            }
                        }
                    )
                else:
                    new_patterns[all_patt_list_name[i][0]]["mols"].append(mol)
                    new_patterns[all_patt_list_name[i][0]]["smiles"].append(smi)
            else:
                max_count = 0
                min_hydrogens = 999
                for j in range(len(all_patt_list[i])):
                    if verbose:
                        print("...... More than one possibility")
                    count = sum(
                        [
                            listElem.count(all_patt_list_name[i][j])
                            for listElem in all_patt_list_name
                        ]
                    )

                    hydrogens = all_patt_list_name[i][j].count("H")
                    if verbose:
                        print(
                            all_patt_list_name[i][j],
                            "occurs",
                            count,
                            "times and features",
                            hydrogens,
                            "hydrogens",
                        )
                    if count > max_count:
                        new_patt = all_patt_list[i][j]
                        new_patt_name = moltosmiles_transition(new_patt, change_dict)
                        max_count = count
                        min_hydrogens = hydrogens
                    elif hydrogens < min_hydrogens:
                        new_patt = all_patt_list[i][j]
                        new_patt_name = moltosmiles_transition(new_patt, change_dict)
                        max_count = count
                        min_hydrogens = hydrogens
                if verbose:
                    print("...... Take pattern", new_patt_name)
                if new_patt_name not in new_patterns:
                    new_patterns.update(
                        {
                            new_patt_name: {
                                "rule": new_patt,
                                "mols": [
                                    mol,
                                ],
                                "smiles": [
                                    smi,
                                ],
                            }
                        }
                    )
                else:
                    new_patterns[new_patt_name]["mols"].append(mol)
                    new_patterns[new_patt_name]["smiles"].append(smi)

    else:
        # Else try to find best combination with most atoms and least hydrogens
        if verbose:
            print("...... Trying to find best combination of extendible atoms:")
        idx_list_all = list(extensions.keys())
        for i in reversed(range(1, len(idx_list_all) + 1)):
            min_hydrogens = 999
            overall_match = False
            subsets = findsubsets(idx_list_all, i)
            for idx_list in subsets:
                if verbose:
                    print("...... Try subsets:", idx_list)
                proposed_new_rule = extend_by_atom(rule, mols[idx], match, idx_list)
                if check_one_extension(proposed_new_rule, max_possible_rules):
                    proposed_new_rule_smiles = moltosmiles_transition(
                        proposed_new_rule, change_dict
                    )
                    overall_match = True
                    hydrogens = proposed_new_rule_smiles.count("H")
                    if hydrogens < min_hydrogens:
                        new_rule = proposed_new_rule
                        new_rule_smiles = proposed_new_rule_smiles
                        min_hydrogens = hydrogens
            if overall_match:
                if verbose:
                    print("...... Accepting pattern", new_rule_smiles)
                new_patterns.update(
                    {
                        new_rule_smiles: {
                            "rule": new_rule,
                            "mols": mols,
                            "smiles": smiles,
                        }
                    }
                )
                break

    return new_patterns


def get_max_possible_rules(possible_extensions, rule, mols, change_dict):
    """
    Takes a dictionary of possible extensions and computes the corresponding enlarged templates as
    RDKit molecule objects.

    Parameters
    ----------
    possible_extensions: dict
        A dictionary of possible extensions for each input molecule, containing the matching indices,
        atom indices of atoms to extend, and their possible extension as string.
    rule: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    mols: List[rdkit.Chem.Mol]
        List of RDKit molecule objects of the molecules or pseudomolecules which are considered in the current branch.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.

    Returns
    -------
    extension_dict: dict
        Dictionary of possibly extensions, containing the atoms indices of the matching atoms, as well as all the
        possible enlarged templates as string and RDKit molecule.

    Examples
    --------
    >>> smis = ["CCOCO", "CCOOC"]
    >>> mols = [Chem.MolFromSmiles(smi) for smi in smis]
    >>> rule = Chem.MolFromSmiles("CCO", sanitize=False)
    >>> possible_extensions = {'CCOCO': {(0, 1, 2): {2: 'OC'}}, 'CCOOC': {(0, 1, 2): {2: 'OO'}}}
    >>> change_dict={"reac": {"atom": {}, "bond": {}}, "prod": {"atom": {}, "bond": {}}}
    >>> get_max_possible_rules(possible_extensions, rule, mols, change_dict)
    {'CCOCO': {(0, 1, 2): {'COCC': <rdkit.Chem.rdchem.Mol at 0x7ffc4d75bee0>}},
     'CCOOC': {(0, 1, 2): {'CCOO': <rdkit.Chem.rdchem.Mol at 0x7ffc4d75b4e0>}}}
    """

    extension_dict = {}
    for i, smi in enumerate(possible_extensions.keys()):
        extension_dict[smi] = {}
        for match in possible_extensions[smi].keys():
            idx_list = list(possible_extensions[smi][match].keys())
            max_rule = extend_by_atom(rule, mols[i], match, idx_list)
            extension_dict[smi][match] = {
                moltosmiles_transition(max_rule, change_dict): max_rule
            }
    return extension_dict


def get_possible_extensions(rule, mols, smiles):
    """
    Searches for all possible extensions of a template taking into account a list of molecules.

    Parameters
    ----------
    rule: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    mols: List[rdkit.Chem.Mol]
        List of RDKit molecule objects of the molecules or pseudomolecules which are considered in the current branch.
    smiles: List[str]
        List of SMILES strings of the molecules or pseudomolecules which are considered in the current branch.

    Returns
    -------
    extension_dict: dict
        A dictionary of possible extensions for each input molecule, containing the matching indices, atom indices
        of atoms to extend, and their possible extension as string.

    Examples
    --------
    For the molecules CCOCO and CCOOC, the rule CCO matches both molecules, and allows for one possible extension at
    the oxygen, adding a carbon for the first molecule, or an oxygen for the second molecule:

    >>> smis = ["CCOCO", "CCOOC"]
    >>> mols = [Chem.MolFromSmiles(smi) for smi in smis]
    >>> rule = Chem.MolFromSmiles("CCO",sanitize=False)
    >>> get_possible_extensions(rule,mols,smis)
    {'CCOCO': {(0, 1, 2): {2: 'OC'}}, 'CCOOC': {(0, 1, 2): {2: 'OO'}}}
    """

    extension_dict = {}
    # For each molecule, get possible extensions with rule atom idx
    for i in range(len(mols)):
        mol = mols[i]
        smi = smiles[i]
        possible_extensions = {}
        matches = mol.GetSubstructMatches(rule)
        for match in matches:
            extension = {}
            for atom in rule.GetAtoms():
                rule_idx = atom.GetIdx()
                rule_neighbors_idx = [
                    neighbor.GetIdx() for neighbor in atom.GetNeighbors()
                ]
                mol_idx = match[rule_idx]
                mol_neighbors_idx = [
                    neighbor.GetIdx()
                    for neighbor in mol.GetAtomWithIdx(mol_idx).GetNeighbors()
                ]
                if len(rule_neighbors_idx) != len(mol_neighbors_idx):
                    # can extend atoms here, find extension pattern:
                    new_rule = extend_by_atom(rule, mol, match, [rule_idx])
                    editable = Chem.EditableMol(new_rule)
                    for idx in range(len(match))[::-1]:
                        if idx != rule_idx:
                            editable.RemoveAtom(idx)
                    pattern = Chem.MolToSmiles(editable.GetMol(), rootedAtAtom=0)
                    extension[rule_idx] = pattern
        possible_extensions[match] = extension
        extension_dict[smi] = possible_extensions
    return extension_dict


def extend_by_atom(patt, m, match, idx_list, ring_extension=False):
    """
    Extends the current rule ('patt') by adding neighbours of the atoms specified in the list of
    atom indices 'idx_list' according to the molecule 'm'.

    Parameters
    ----------
    patt: rdkit.Chem.Mol
        RDKit molecule object of the current template (molecule or pseudomolecule).
    n: rdkit.Chem.Mol
        RDKit molecule object of the current leaf node.
    match: tuple
        Tuple of matching atom indices of the current rule.
    idx_list: list
        List of atom indices at which to extend the pattern.
    ring_extension: bool
        Boolean whether current iteration has broken a ring and must thus be iterated until full ring is found.

    Returns
    -------
    new_patt: rdkit.Chem.Mol
        RDKit molecule object of the new, extended template (molecule or pseudomolecule).

    Examples
    --------
    For the molecule CCOCO, the rule CCO matches the first three atoms and can be extended at atom 2 (the oxygen),
    yielding the new, extended molecule CCOC:

    >>> rule = Chem.MolFromSmiles("CCO",sanitize=False)
    >>> new_rule=extend_by_atom(rule,Chem.MolFromSmiles("CCOCO"),(0, 1, 2),[2])
    >>> print(Chem.MolToSmiles(new_rule))
    'COCC'
    """

    added_atoms = {}
    ring = False
    new_patt = deepcopy(patt)
    for pattern_idx in idx_list:
        match_idx = match[pattern_idx]
        for neighbor in m.GetAtomWithIdx(match_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if (
                neighbor_idx not in match
            ):  # Need to add new bond (and possibly new atom if not already added)
                if neighbor_idx not in added_atoms.keys():  # New bond and atom
                    if m.GetBondBetweenAtoms(match_idx, neighbor_idx).IsInRing():
                        ring = True
                    if (
                        ring_extension
                        and m.GetBondBetweenAtoms(match_idx, neighbor_idx).IsInRing()
                        or not ring_extension
                    ):  # If in ring extension mode, only add ring atoms
                        editable = Chem.EditableMol(new_patt)
                        idx = editable.AddAtom(Chem.Atom(neighbor.GetAtomicNum()))
                        editable.AddBond(
                            pattern_idx,
                            idx,
                            order=m.GetBondBetweenAtoms(
                                match_idx, neighbor_idx
                            ).GetBondType(),
                        )
                        new_patt = editable.GetMol()
                        new_patt.GetAtomWithIdx(idx).SetNoImplicit(True)
                        new_patt.GetAtomWithIdx(idx).SetFormalCharge(
                            m.GetAtomWithIdx(neighbor_idx).GetFormalCharge()
                        )
                        added_atoms[neighbor_idx] = idx
                        match = match + (neighbor_idx,)
                else:  # found a ring, need to reconnect it, add only new bond (not atom)
                    editable = Chem.EditableMol(new_patt)
                    idx = added_atoms[neighbor_idx]
                    editable.AddBond(
                        pattern_idx,
                        idx,
                        order=m.GetBondBetweenAtoms(
                            match_idx, neighbor_idx
                        ).GetBondType(),
                    )
                    new_patt = editable.GetMol()
            elif not new_patt.GetBondBetweenAtoms(
                pattern_idx, match.index(neighbor_idx)
            ):  # bond between known atoms missing
                editable = Chem.EditableMol(new_patt)
                editable.AddBond(
                    pattern_idx,
                    match.index(neighbor_idx),
                    order=m.GetBondBetweenAtoms(match_idx, neighbor_idx).GetBondType(),
                )
                new_patt = editable.GetMol()
    if ring:
        new_patt = extend_by_atom(
            new_patt, m, match, list(added_atoms.values()), ring_extension=True
        )

    return new_patt
