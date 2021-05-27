"""
transition_state.py
Some of the functions in this file originate from RDChiral (https://github.com/connorcoley/rdchiral).
"""

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from copy import deepcopy
import re
import subprocess
import os


def make_mols(smi, stereoinformation):
    """
    Takes a smiles string as input and creates an RDKit mol object. The molecule cannot be sanitized,
    since we would else get rid of the explicit hydrogens. Thus, instead we check if the molecule could
    be sanitized, and reorder the atoms corresponding to the sanitized molecule. The reordering of atoms
    and molecules ensures that the same reaction with different atom ordering gives the same transition
    state.

    Parameters
    ----------
    smi: str
        SMILES string.
    stereoinformation: bool
        Whether to use stereoinformation.

    Returns
    -------
    mol_ordered: rdkit.Chem.Mol
        RDKit molecule object.
    """

    if Chem.MolFromSmiles(smi):
        # Order molecules
        sorted_smi = ""
        smis = smi.split(".")
        smis_nolabel = re.sub(r"\:[0-9]+\]", "]", smi).split(".")
        order = []
        natoms = []
        for i in range(len(smis_nolabel)):
            if not stereoinformation:
                smis[i] = (
                    smis[i].replace("/", "").replace("\\", "").replace("@", "")
                )  # Delete Stereo
            mol = Chem.MolFromSmiles(smis[i], sanitize=False)
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
            )  # Sanitize all functions apart hydrogens
            mol_nolabel = Chem.MolFromSmiles(smis_nolabel[i])
            smis_nolabel[i] = Chem.MolToSmiles(mol_nolabel)
            order.append(
                list(
                    mol.GetSubstructMatch(
                        Chem.AddHs(
                            Chem.MolFromSmiles(
                                Chem.MolToSmiles(
                                    Chem.AddHs(Chem.MolFromSmiles(smis_nolabel[i]))
                                )
                            )
                        )
                    )
                )
            )
            natoms.append(mol.GetNumAtoms())
        overall_order = []
        sorting_keys = sorted(range(len(smis_nolabel)), key=lambda k: smis_nolabel[k])
        j = 0
        for key in sorting_keys:
            sorted_smi += smis[key]
            sorted_smi += "."
            overall_order.extend([x + j for x in order[key]])
            j += natoms[key]
        sorted_smi = sorted_smi[:-1]

        # Create molecule and reorder atom
        mol = Chem.MolFromSmiles(sorted_smi, sanitize=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
        )  # Sanitize all functions apart hydrogens
        mol_ordered = Chem.RenumberAtoms(mol, overall_order)
        mol_ordered = Chem.MolFromSmiles(
            Chem.MolToSmiles(mol_ordered), sanitize=False
        )  # Correct atom ordering needs another round of mol -> smiles -> mol
        Chem.SanitizeMol(
            mol_ordered,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
        )  # Sanitize all functions apart hydrogens
        mol_ordered.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol_ordered)
    else:
        raise ValueError("Smiles string doesn't produce a valid RDKit mol!")

    return mol_ordered


def bond_to_label(bond):
    """
    This function takes an RDKit bond and creates a label describing
    the most important attributes.

    Parameters
    ----------
    bond: rdkit.Chem.Bond
        RDKit bond.

    Returns
    -------
    label: str
        Label of bond.
    """

    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().HasProp("molAtomMapNumber"):
        a1_label += bond.GetBeginAtom().GetProp("molAtomMapNumber")
    if bond.GetEndAtom().HasProp("molAtomMapNumber"):
        a2_label += bond.GetEndAtom().GetProp("molAtomMapNumber")
    atoms = sorted([a1_label, a2_label])

    label = "{}{}{}".format(atoms[0], bond.GetSmarts(), atoms[1])

    return label


def atoms_are_different(atom1, atom2):
    """
    Compares two RDKit atoms based on basic properties.

    Parameters
    ----------
    atom1: rdkit.Chem.Atom
        First RDKit atom.
    atom2: rdkit.Chem.Atom
        Second RDKit atom.

    Returns
    -------
    bool
        Boolean whether the atoms have the same properties.
    """

    if atom1.GetSmarts() != atom2.GetSmarts():
        return True  # should be very general
    if atom1.GetAtomicNum() != atom2.GetAtomicNum():
        return True  # must be true for atom mapping
    if atom1.GetFormalCharge() != atom2.GetFormalCharge():
        return True
    if atom1.GetDegree() != atom2.GetDegree():
        return True
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons():
        return True
    if atom1.GetIsAromatic() != atom2.GetIsAromatic():
        return True

    # Check bonds and nearest neighbor identity
    bonds1 = sorted([bond_to_label(bond) for bond in atom1.GetBonds()])
    bonds2 = sorted([bond_to_label(bond) for bond in atom2.GetBonds()])
    if bonds1 != bonds2:
        return True

    return False


def get_tagged_atoms_from_mol(mol):
    """Takes an RDKit molecule and returns list of tagged atoms and their
    corresponding numbers.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.

    Returns
    -------
    atoms: List[rdkit.Chem.Atom]
        List of tagged atoms
    atom_tags: List[str]
        List of atom-mapping numbers
    """

    atoms = []
    atom_tags = []
    for atom in mol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            atoms.append(atom)
            atom_tags.append(str(atom.GetProp("molAtomMapNumber")))
    return atoms, atom_tags


def get_changed_atoms(reacs, prods):
    """
    Looks at mapped atoms in a reaction and determines which ones changed.

    Parameters
    ----------
    reacs: rdkit.Chem.Mol
        RDKit molecule of reactant(s).
    prods: rdkit.Chem.Mol
        RDKit molecule of product(s).

    Returns
    -------
    changed_atoms: List[rdkit.Chem.Atom]
        List of changed atoms
    changed_atom_tags: List[str]
        List of tag numbers of changed atoms
    err: int
        Integer indicating an error if not equal 0.
    """

    err = 0
    prod_atoms, prod_atom_tags = get_tagged_atoms_from_mol(prods)
    reac_atoms, reac_atom_tags = get_tagged_atoms_from_mol(reacs)
    if len(set(prod_atom_tags)) != len(set(reac_atom_tags)):
        print("warning: different atom tags appear in reactants and products")
        err = 1
    if len(prod_atoms) != len(reac_atoms):
        print("warning: total number of tagged atoms differ, stoichometry != 1?")
        err = 1

    # Find differences
    changed_atoms = []
    changed_atom_tags = []

    # Product atoms that are different from reactant atom equivalent
    for i, prod_tag in enumerate(prod_atom_tags):
        for j, reac_tag in enumerate(reac_atom_tags):
            if reac_tag != prod_tag:
                continue
            if (
                reac_tag not in changed_atom_tags
            ):  # don't bother comparing if we know this atom changes
                # If atom changed, add
                if atoms_are_different(prod_atoms[i], reac_atoms[j]):
                    changed_atoms.append(reac_atoms[j])
                    changed_atom_tags.append(reac_tag)
                    break

    return changed_atoms, changed_atom_tags, err


def get_strict_smarts_for_atom(atom):
    """
    For an RDkit atom object, generate a simple SMARTS pattern.

    Parameters
    ----------
    atom: rdkit.Chem.Atom
        RDKit atom.

    Results
    -------
    symbol: str
        SMARTS symbol of the atom.
    """

    symbol = atom.GetSmarts()
    if "[" not in symbol:
        symbol = "[" + symbol + "]"

    return symbol


def get_fragments_for_changed_atoms(mol, changed_atom_tags):
    """
    Given an RDKit mols and a list of changed atom tags, this function
    computes the SMILES string of molecular fragments using MolFragmentToSmiles
    for all changed fragments.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule
    changed_atom_tags: List[str]
        List of changed atom tags.

    Returns
    -------
    this_fragment: str
        SMILES string of molecular fragment of changed atoms
    fragment_atoms: dict
        Dictionary of atoms in fragment.
    """

    fragment_atoms = {}

    # Initialize list of replacement symbols (updated during expansion)
    symbol_replacements = []

    # Build list of atoms to use
    atoms_to_use = []
    for atom in mol.GetAtoms():
        # Check self (only tagged atoms)
        if ":" in atom.GetSmarts():
            if atom.GetSmarts().split(":")[1][:-1] in changed_atom_tags:
                fragment_atoms[int(atom.GetSmarts().split(":")[1][:-1])] = atom.GetIdx()
                atoms_to_use.append(atom.GetIdx())
                symbol = get_strict_smarts_for_atom(atom)
                if symbol != atom.GetSmarts():
                    symbol_replacements.append((atom.GetIdx(), symbol))
                continue

    # Define new symbols based on symbol_replacements
    symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
    for (i, symbol) in symbol_replacements:
        symbols[i] = symbol

    mol_copy = deepcopy(mol)
    [x.ClearProp("molAtomMapNumber") for x in mol_copy.GetAtoms()]
    this_fragment = AllChem.MolFragmentToSmiles(
        mol_copy,
        atoms_to_use,
        atomSymbols=symbols,
        allHsExplicit=True,
        isomericSmiles=False,
        allBondsExplicit=True,
    )
    return this_fragment, fragment_atoms


def changed_labels(
    reac_fragments,
    prod_fragments,
    reac_fragment_atoms,
    prod_fragment_atoms,
    reac_mols,
    prod_mols,
):
    """
    Given the fragments of reactants and products (strings), the mapping of index to atom-mapping number
    (dictionaries), and the whole molecules (RDKit mols), this function computes the changed labels, this
    is the changed formal charges (e.g. R-O-H to R-O- (change of charge on the oxygen).

    Parameters
    ----------
    reac_fragments: str
        Fragment of reactant(s).
    prod_fragments: str
        Fragment of product(s).
    reac_fragment_atoms: dict
        Dictionary of index to atom-mapping numbers in reactant fragment.
    prod_fragment_atoms: dict
        Dictionary of index to atom-mapping numbers in product fragment.
    reac_mols: rdkit.Chem.Mol
        RDKit reactant(s) molecule.
    reac_prod: rdkit.Chem.Mol
        RDKit product(s) molecule.

    Returns
    change_dict_atoms: list
        List of changed labels (=changed formal charges).
    """

    change_dict_atoms = {}
    tags = re.findall(r"\:([0-9]+)\]", reac_fragments)
    reac_symbols = re.findall(r"\[[^\]]*\]", reac_fragments)
    prod_symbols = re.findall(r"\[[^\]]*\]", prod_fragments)

    for tag in tags:
        reac = reac_symbols[
            [(":" + tag + "]" in reac_symbol) for reac_symbol in reac_symbols].index(
                True
            )
        ]
        prod = prod_symbols[
            [(":" + tag + "]" in prod_symbol) for prod_symbol in prod_symbols].index(
                True
            )
        ]
        if reac != prod:
            charge_reac = reac_mols.GetAtomWithIdx(
                reac_fragment_atoms[int(tag)]
            ).GetFormalCharge()
            charge_prod = prod_mols.GetAtomWithIdx(
                prod_fragment_atoms[int(tag)]
            ).GetFormalCharge()
            change_dict_atoms[int(tag)] = [charge_reac, charge_prod]

    return change_dict_atoms


def include_bonds(mols, changed_atom_tags, fragment_atoms):
    """Creates a list of bonds and their indices that have changed atoms attached and should thus
    be included in the transition state.

    Parameters
    ----------
    mols: List[rdkit.Chem.Mol]
        RDKit molecule(s).
    changed_atom_tags: List[int]
        Indices of changed atoms.
    fragment_atoms: dict
        Dictionary of index to atom-mapping numbers in fragment.

    Returns
    -------
    included_bonds: list
        List of included bonds.
    included_bonds_idx: List[int]
        List of included bond indices
    """

    included_bonds = []
    included_bonds_idx = []
    for idx in fragment_atoms.values():
        atom = mols.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (
                begin_atom.GetProp("molAtomMapNumber") in changed_atom_tags
                and end_atom.GetProp("molAtomMapNumber") in changed_atom_tags
            ):
                if (
                    idx <= begin_atom.GetIdx() and idx <= end_atom.GetIdx()
                ):  # cound each bond only once
                    included_bonds.append(
                        [
                            int(begin_atom.GetProp("molAtomMapNumber")),
                            int(end_atom.GetProp("molAtomMapNumber")),
                            bond.GetBondType(),
                        ]
                    )
                    included_bonds_idx.append(bond.GetIdx())
    return included_bonds, included_bonds_idx


def changed_bonds(included_bonds_reac, included_bonds_prod):
    """
    Creates a dictionary of changed bonds from a list of included bonds in the reactants and products.

    Parameters
    ----------

    included_bonds_reac: list
        List of included bonds in the reactants.
    included_bonds_prod: list
        List of included bonds in the products.

    Returns
    -------
    change_dict_bonds: dict
        Dictionary of changed bonds.
    """

    change_dict_bonds = {}
    for bond in included_bonds_reac:
        if bond not in included_bonds_prod:
            change_dict_bonds[tuple(sorted([bond[0], bond[1]]))] = [
                bond[2],
                Chem.rdchem.BondType.ZERO,
            ]
    for bond in included_bonds_prod:
        if bond not in included_bonds_reac:
            if tuple(sorted([bond[0], bond[1]])) in change_dict_bonds.keys():
                change_dict_bonds[tuple(sorted([bond[0], bond[1]]))][1] = bond[2]
            else:
                change_dict_bonds[tuple(sorted([bond[0], bond[1]]))] = [
                    Chem.rdchem.BondType.ZERO,
                    bond[2],
                ]

    return change_dict_bonds


def bond_order_increases(bond1, bond2):
    """
    Determines whether bond order gets larger from bond1 to bond2.

    Parameters
    ----------
    bond1: rdkit.Chem.Bond
        RDKit bond object.
    bond2: rdkit.Chem.Bond
        RDKit bond object.

    Returns
    -------
    bool
        Whether bond order increases from bond1 to bond2.
    """

    known_bonds = [
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DOUBLE,
        Chem.rdchem.BondType.TRIPLE,
        Chem.rdchem.BondType.AROMATIC,
        Chem.rdchem.BondType.ZERO,
    ]
    assert (
        bond1 in known_bonds
    ), "Bond type not supported (must be zero, single, double, triple or aromatic)"
    assert (
        bond2 in known_bonds
    ), "Bond type not supported (must be zero, single, double, triple or aromatic)"

    values = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.AROMATIC: 1.5,
        Chem.rdchem.BondType.ZERO: 0,
    }
    if values[bond2] > values[bond1]:
        return True
    else:
        return False


def make_ts(
    reac_mols,
    included_bonds_idx,
    change_dict_bonds,
    change_dict_atoms,
    reac_fragment_atoms,
    ts_seed,
):
    """
    Makes an RDKit mol of the transition state, either the whole transistion state (with all atoms and bonds included,
    or only the seed (if ts_seed==True), which includes only the atoms undergoing changes and the corresponding bonds.
    Returns the transition state (whole or seed), and dictionaries of how to change bond types and charges to
    revert back to either the products or reactants.

    Parameters
    ----------

    reac_mols: rdkit.Chem.mol
        RDKit molecule of reactant(s).
    included_bonds_idx: List[int]
        List of included bond indices in transition state.
    change_dict_bonds: dict
        Dictionary of bond changes in reaction.
    change_dict_atoms: dict
        Dictionary of formal charge changes in reaction.
    reac_fragment_atoms: dict
        Dictionary of atom indices in reactant(s) fragment.
    ts_seed: bool
        Boolean whether only to save the full transition state or only the reaction center.

    Returns
    -------
    change_ts_to_reac: dict
        Dictionary of bond and atom changes from the ts to the reactants
    change_ts_to_prod: dict
        Dictionary of bond and atom changes from the ts to the products,
    ts: rdkit.Chem.Mol
        RDKit molecule of the transition state
    """

    ts = deepcopy(reac_mols)
    ts_included_bonds_idx = deepcopy(included_bonds_idx)
    bond_idx = ts.GetNumBonds()
    for key in change_dict_bonds.keys():
        if (
            change_dict_bonds[key][0] != Chem.rdchem.BondType.ZERO
        ):  # Change of an existing bond
            bond = ts.GetBondBetweenAtoms(
                reac_fragment_atoms[key[0]], reac_fragment_atoms[key[1]]
            )
            if bond_order_increases(
                change_dict_bonds[key][0], change_dict_bonds[key][1]
            ):
                bond.SetBondType(Chem.rdchem.BondType.DATIVEL)
            else:
                bond.SetBondType(Chem.rdchem.BondType.DATIVER)
        else:  # add new bond
            editable = Chem.EditableMol(ts)
            editable.AddBond(
                reac_fragment_atoms[key[0]],
                reac_fragment_atoms[key[1]],
                order=Chem.rdchem.BondType.DATIVEL,
            )
            ts_included_bonds_idx.append(bond_idx)
            ts = editable.GetMol()
            bond_idx += 1

    if ts_seed:
        amap = {}
        ts = Chem.PathToSubmol(
            ts, tuple(ts_included_bonds_idx), atomMap=amap
        )  # Match only fragment
    else:
        amap = {}
        ts = Chem.PathToSubmol(
            ts, tuple(range(bond_idx)), atomMap=amap
        )  # Match whole molecule

    # Set up dictionary of change
    change_ts_to_reac = {"bond": {}, "atom": {}}
    change_ts_to_prod = {"bond": {}, "atom": {}}
    for key in change_dict_bonds.keys():
        id_atom1 = amap[reac_fragment_atoms[key[0]]]
        id_atom2 = amap[reac_fragment_atoms[key[1]]]
        bond_idx = (id_atom1, id_atom2)
        if bond_order_increases(change_dict_bonds[key][0], change_dict_bonds[key][1]):
            change_ts_to_reac["bond"][bond_idx] = (
                change_dict_bonds[key][0],
                Chem.rdchem.BondType.DATIVEL,
            )
            change_ts_to_prod["bond"][bond_idx] = (
                change_dict_bonds[key][1],
                Chem.rdchem.BondType.DATIVEL,
            )
        else:
            change_ts_to_reac["bond"][bond_idx] = (
                change_dict_bonds[key][0],
                Chem.rdchem.BondType.DATIVER,
            )
            change_ts_to_prod["bond"][bond_idx] = (
                change_dict_bonds[key][1],
                Chem.rdchem.BondType.DATIVER,
            )

    for key in change_dict_atoms.keys():
        id_atom = amap[reac_fragment_atoms[key]]
        change_ts_to_reac["atom"][id_atom] = (
            change_dict_atoms[key][0],
            change_dict_atoms[key][0],
        )
        change_ts_to_prod["atom"][id_atom] = (
            change_dict_atoms[key][1],
            change_dict_atoms[key][0],
        )

    return change_ts_to_reac, change_ts_to_prod, ts


def mapping_list_dict(mol):
    """
    Computes a list of all map numbers and an atom-mapping indices dictionary.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.

    Returns
    -------
    map_numbers: List[int]
        Sorted list of atom map numbers
    mapping: dict
        Dictionary of mapping
    """

    map_numbers = []
    mapping = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            map_numbers.append(int(atom.GetProp("molAtomMapNumber")))
            mapping[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    map_numbers.sort()

    return map_numbers, mapping


def compute_aam_with_h(rxn_smiles):
    """
    Helper function to produce an atom mapping for a reaction smiles via the
    ReactionRecoder tool (RDT) by Syed Asad Rahman. The tool maps hydrogens that undergo changes,
    and this function adds the necessary atom-mapping to all other hydrogens.

    Parameters
    ----------
    rxn_smiles: str
        Reaction SMILES.

    Returns
    -------
    rxn_smiles: str
        atom-mapped reaction SMILES with hydrogens.
    """

    reac_smiles = Chem.MolToSmiles(
        Chem.AddHs(Chem.MolFromSmiles(rxn_smiles.split(">")[0]))
    )
    prod_smiles = Chem.MolToSmiles(
        Chem.AddHs(Chem.MolFromSmiles(rxn_smiles.split(">")[-1]))
    )
    rxn_smiles = reac_smiles + ">>" + prod_smiles

    dirname = os.path.dirname(__file__)
    rdt_path = os.path.join(dirname, "../../ReactionDecoder/ReactionDecoder.jar")

    s = f'java -jar {rdt_path} -Q SMI -q "{rxn_smiles}"  -j AAM -f TEXT'
    subprocess.run(s, shell=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        'grep -A 1 "SELECTED AAM MAPPING" ECBLAST_smiles_AAM.txt | tail -n 1 > ECBLAST_smiles_AAM.smi',
        shell=True,
    )
    with open("ECBLAST_smiles_AAM.smi") as f:
        rxn_smiles = f.read().splitlines()
    subprocess.check_call("rm -f ECBLAST_smiles_AAM.*", shell=True)
    reac_mols = Chem.AddHs(
        Chem.MolFromSmiles(rxn_smiles[0].split(">")[0], sanitize=False)
    )
    prod_mols = Chem.AddHs(
        Chem.MolFromSmiles(rxn_smiles[0].split(">")[-1], sanitize=False)
    )
    Chem.SanitizeMol(
        reac_mols,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
    )  # Sanitize all functions apart hydrogens
    Chem.SanitizeMol(
        prod_mols,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
    )  # Sanitize all functions apart hydrogens

    map_numbers_reac, mapping_reac = mapping_list_dict(reac_mols)
    map_numbers_prod, mapping_prod = mapping_list_dict(prod_mols)

    if map_numbers_reac != map_numbers_prod:
        raise ValueError("Atom mapping is missing atoms. Exit program.")

    idx = int(np.max(map_numbers_reac)) + 1
    for num in map_numbers_reac:
        atom_reac = reac_mols.GetAtomWithIdx(mapping_reac[num])
        atom_prod = prod_mols.GetAtomWithIdx(mapping_prod[num])
        if atom_reac.GetSymbol() != "H":
            idx2 = 0
            for atom_reac_neighbor in atom_reac.GetNeighbors():
                if (
                    atom_reac_neighbor.GetSymbol() == "H"
                    and not atom_reac_neighbor.HasProp("molAtomMapNumber")
                ):
                    atom_reac_neighbor.SetAtomMapNum(idx + idx2)
                    idx2 += 1
            idx3 = 0
            for atom_prod_neighbor in atom_prod.GetNeighbors():
                if (
                    atom_prod_neighbor.GetSymbol() == "H"
                    and not atom_prod_neighbor.HasProp("molAtomMapNumber")
                ):
                    atom_prod_neighbor.SetAtomMapNum(idx + idx3)
                    idx3 += 1
            if idx2 != idx3:
                raise ValueError("Atom mapping is missing atoms. Exit program.")
            idx += idx2

    rxn_smiles = Chem.MolToSmiles(reac_mols) + ">>" + Chem.MolToSmiles(prod_mols)

    return rxn_smiles


def compute_aam_without_h(rxn_smiles):
    """
    This helper function produces an atom mapping for a reaction smiles via the
    ReactionRecoder tool (RDT) by Syed Asad Rahman.

    Parameters
    ----------
    rxn_smiles: str
        Reaction SMILES.

    Returns
    -------
    rxn_smiles: str
        atom-mapped reaction SMILES without hydrogens.
    """

    reac_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(rxn_smiles.split(">")[0]))
    prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(rxn_smiles.split(">")[-1]))
    rxn_smiles = reac_smiles + ">>" + prod_smiles

    dirname = os.path.dirname(__file__)
    rdt_path = os.path.join(dirname, "../../ReactionDecoder/ReactionDecoder.jar")

    s = f'java -jar {rdt_path} -Q SMI -q "{rxn_smiles}"  -j AAM -f TEXT'
    subprocess.run(s, shell=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        'grep -A 1 "SELECTED AAM MAPPING" ECBLAST_smiles_AAM.txt | tail -n 1 > ECBLAST_smiles_AAM.smi',
        shell=True,
    )
    with open("ECBLAST_smiles_AAM.smi") as f:
        rxn_smiles = f.read().splitlines()
    subprocess.check_call("rm -f ECBLAST_smiles_AAM.*", shell=True)

    return rxn_smiles[0]


def process_an_example(rxn_smiles):
    """
    Processes one rxn_smiles string and returns the transition state (whole or seed depending on ts_seed)
    RDKit mol object, as well as dictionaries of how to change bond types and charges to revert the transition state
    back to either the products or reactants.

    Parameters
    ----------
    rxn_smiles: str
        Reaction SMILES.

    Returns
    -------
    change_ts_to_reac: dict
        Dictionary of bond and atom changes from the ts to the reactants
    change_ts_to_prod: dict
        Dictionary of bond and atom changes from the ts to the products
    ts: rdkit.Chem.Mol
        RDKit molecule of the whole imaginary transition state
    ts_seed: rdkit.Chem.Mol
        RDKit molecule of the reaction center
    """

    reac_mols = make_mols(rxn_smiles.split(">")[0], False)
    prod_mols = make_mols(rxn_smiles.split(">")[-1], False)

    changed_atoms, changed_atom_tags, err = get_changed_atoms(reac_mols, prod_mols)
    if err:
        raise ValueError("Could not get changed atoms")
    if not changed_atom_tags:
        raise ValueError("No change detected")

    reac_fragments, reac_fragment_atoms = get_fragments_for_changed_atoms(
        reac_mols, changed_atom_tags
    )
    prod_fragments, prod_fragment_atoms = get_fragments_for_changed_atoms(
        prod_mols, changed_atom_tags
    )
    included_bonds_reac, included_bonds_idx = include_bonds(
        reac_mols, changed_atom_tags, reac_fragment_atoms
    )
    included_bonds_prod, _ = include_bonds(
        prod_mols, changed_atom_tags, prod_fragment_atoms
    )
    change_dict_bonds = changed_bonds(included_bonds_reac, included_bonds_prod)
    change_dict_atoms = changed_labels(
        reac_fragments,
        prod_fragments,
        reac_fragment_atoms,
        prod_fragment_atoms,
        reac_mols,
        prod_mols,
    )

    change_ts_to_reac, change_ts_to_prod, ts_seed = make_ts(
        reac_mols,
        included_bonds_idx,
        change_dict_bonds,
        change_dict_atoms,
        reac_fragment_atoms,
        True,
    )
    _, _, ts = make_ts(
        reac_mols,
        included_bonds_idx,
        change_dict_bonds,
        change_dict_atoms,
        reac_fragment_atoms,
        False,
    )

    return change_ts_to_reac, change_ts_to_prod, ts, ts_seed
