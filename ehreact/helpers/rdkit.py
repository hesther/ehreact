from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from copy import deepcopy
import re
import itertools
from ehreact.helpers.transition_state import process_an_example


def morgan_bit_fp_from_mol(mol, chirality=False):
    """
    Computes a Morgan Bit Fingerprint for a molecule.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    chirality: bool
        Whether to use stereoinformation in fingerprint.

    Returns
    -------
    fp: rdkit.DataStructs.cDataStructs.ExplicitBitVect
        Morgan fingerprint.
    """

    mol = Chem.RemoveHs(mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=2048, useFeatures=False, useChirality=chirality
    )
    return fp


def tanimoto_from_fp(fp1, fp2):
    """
    Computes the Tanimoto similarity between fingerprints.

    Parameters
    ----------
    fp1: rdkit.DataStructs.cDataStructs.ExplicitBitVect
        Molecular fingerprint.
    fp2: rdkit.DataStructs.cDataStructs.ExplicitBitVect
        Molecular fingerprint.

    Returns
    -------
    tanimoto: float
        Tanimoto similarity
    """

    tanimoto = DataStructs.cDataStructs.TanimotoSimilarity(fp1, fp2)
    return tanimoto


def find_common(mols):
    """
    Function to find the maximum common substructure of a list of RDKit molecules
    and check whether charges match.

    Parameters
    ----------
    mols: List[rdkit.Chem.Mol]
        List of RDKit mols

    Returns
    -------
    new_mcs_mol: rdkit.Chem.Mol
        RDKit molecule of common substructure.
    """

    if len(mols) <= 1:
        raise ValueError(
            "Not enough molecules to calculate common structures. Provide a seed instead."
        )

    atomComp = rdFMCS.AtomCompare.CompareElements

    atomParam = rdFMCS.MCSAtomCompareParameters()
    atomParam.MatchChiralTag = False
    atomParam.MatchFormalCharge = True
    atomParam.MatchValences = False
    atomParam.RingMatchesRingOnly = False

    bondComp = rdFMCS.BondCompare.CompareOrderExact

    bondParam = rdFMCS.MCSBondCompareParameters()
    bondParam.CompleteRingsOnly = True
    bondParam.RingMatchesRingOnly = False

    opt = rdFMCS.MCSParameters()
    opt.MaximizeBonds = True
    opt.Timeout = 3600
    opt.AtomCompareParameters = atomParam
    opt.BondCompareParameters = bondParam
    opt.SetAtomTyper(atomComp)
    opt.SetBondTyper(bondComp)
    opt.Threshold = 1.0
    opt.Verbose = False
    opt.seedSmarts = ""

    mcs = rdFMCS.FindMCS(mols, opt)

    if mcs.canceled:
        print("WARNING: MCS search timed out, returning best match")
    if mcs.numAtoms == 0:
        raise ValueError("No MCS was found between the molecules")

    # Convert to RDKit molecule:
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # Reintroduce charges:
    matches = mols[0].GetSubstructMatches(mcs_mol)
    new_mcs_mol = None
    for match in matches:
        mcs_mol_charged = force_charge_fit(mols[0], mcs_mol, match)

        # check if charged mcs mol is ok with other mols
        ok = []
        for mol in mols[1:]:
            matches2 = mol.GetSubstructMatches(mcs_mol_charged)
            if True in [
                do_charges_fit(mol, mcs_mol_charged, match2) for match2 in matches2
            ]:
                ok.append(True)
                break
            ok.append(False)
        if False not in ok:
            new_mcs_mol = mcs_mol_charged
            break
    if new_mcs_mol is None:
        raise ValueError("Could not find consistent charges for MCS")

    return new_mcs_mol


def check_if_reaction(smi, expected):
    """
    Function to check whether ">" occurs in a string, and raise a ValueError if the expected behavior is not met.

    Parameters
    ----------
    smi: str
        SMILES string (reaction SMILES or molecule SMILES)
    expected: bool
        Whether the input string is expected to be a reaction.
    """

    observed = ">" in smi
    if observed and not expected:
        raise ValueError("Expected molecule smiles but received reaction smiles.")
    elif expected and not observed:
        raise ValueError("Expected reaction smiles but received molecule smiles.")


def make_mol(smi):
    """
    Make RDKit molecule from a smiles string and add hydrogens.

    Parameters
    ----------
    smi: str
        SMILES string.

    Returns
    -------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    """

    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    return mol


def make_mol_no_sanit_h(smi):
    """
    Make RDKit molecule from a smiles string without sanitizing hydrogens (keep hydrogens as given).

    Parameters
    ----------
    smi: str
        SMILES string.

    Returns
    -------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    """

    mol = Chem.MolFromSmiles(smi, sanitize=False)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
    )  # Sanitize all functions apart hydrogens
    return mol


def canonicalize(mol):
    """
    Outputs the canonical SMILES string of a molecule.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.

    Returns
    -------
    canonical_smi: str
        Canonical SMILES string.
    """

    smi = Chem.MolToSmiles(mol, canonical=True)
    mol = make_mol(smi)
    canonical_smi = Chem.MolToSmiles(mol, canonical=True)
    return canonical_smi


def canonicalize_with_h(mol, reaction=False, old2new_mapno={}):
    """
    Outputs the canonical SMILES with hydrogens.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    reaction: bool, default False
        Whether the input SMILES stems from a reaction and therefore map numbers need to be included
    old2new_mapno: dict, default {}
        Dictionary of custom map numbers to use. If not supplied, map numbers will be calculated from
        atom indices.

    Returns
    -------
    canonical_smi: str
        Canonical SMILES string.
    old2new_mapno: dict
        Dictionary of map numbers.
    """
    if reaction:
        if old2new_mapno == {}:
            for atom in mol.GetAtoms():
                old2new_mapno[atom.GetAtomMapNum()] = atom.GetIdx() + 1
                atom.SetAtomMapNum(atom.GetIdx() + 1)
        else:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(old2new_mapno[atom.GetAtomMapNum()])

    smi = Chem.MolToSmiles(mol, canonical=True)
    mol = make_mol_no_sanit_h(smi)
    canonical_smi = Chem.MolToSmiles(mol, canonical=True)

    return canonical_smi, old2new_mapno


def preprocess(smi, delete_aam=True, reaction=False, old2new_mapno={}):
    """
    Preprocess a SMILES string.

    Parameters
    ----------
    smi: str
        SMILES string.
    delete_aam: bool, default True
        Whether to delete atom mappings from the string.
    reaction: bool, default False
        Whether the input SMILES stems from a reaction and molecules need to be created
        without sanitizing hydrogens.
    old2new_mapno: dict, default {}
        Dictionary of custom map numbers to use.

    Returns
    -------
    canonical_smi_no_stereo: str
        Canonical SMILES string without stereoinformation.
    canonical_smi: str
        Canonical SMILES string with stereoinformation (if provided).
    mol_no_stereo: rdkit.Chem.Mol
        RDKit molecule without stereoinformation.
    mol: rdkit.Chem.Mol
        RDKit molecule with stereoinformation.
    old2new_mapno: dict
        Dictionary of map numbers.
    """

    if delete_aam:
        smi = re.sub(r"\:[0-9]+\]", "]", smi)

    # Without stereochemistry:
    smi_no_stereo = smi.replace("/", "").replace("\\", "").replace("@", "")
    if reaction:
        mol_no_stereo = make_mol_no_sanit_h(smi_no_stereo)
    else:
        mol_no_stereo = Chem.AddHs(make_mol(smi_no_stereo))
    canonical_smi_no_stereo, _ = canonicalize_with_h(
        mol_no_stereo, reaction, old2new_mapno
    )

    # With stereochemistry:
    if reaction:
        mol = make_mol_no_sanit_h(smi)
    else:
        mol = Chem.AddHs(make_mol(smi))
    canonical_smi, old2new_mapno = canonicalize_with_h(mol, reaction, old2new_mapno)

    return canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol, old2new_mapno


def preprocess_rxn(rxn_smi):
    """
    Preprocess a reaction SMILES string.

    Parameters
    ----------
    rxn_smi: str
        Reaction SMILES string.

    Returns
    -------
    canonical_smi_no_stereo: str
        Canonical SMILES string without stereoinformation.
    canonical_smi: str
        Canonical SMILES string with stereoinformation (if provided).
    mol_no_stereo: rdkit.Chem.Mol
        RDKit molecule without stereoinformation.
    mol: rdkit.Chem.Mol
        RDKit molecule with stereoinformation.
    """

    reac = rxn_smi.split(">")[0]
    (
        reac_canonical_smi_no_stereo,
        reac_canonical_smi,
        reac_mol_no_stereo,
        reac_mol,
        old2new_mapno,
    ) = preprocess(reac, delete_aam=False, reaction=True, old2new_mapno={})

    prod = rxn_smi.split(">")[-1]
    (
        prod_canonical_smi_no_stereo,
        prod_canonical_smi,
        prod_mol_no_stereo,
        prod_mol,
        _,
    ) = preprocess(prod, delete_aam=False, reaction=True, old2new_mapno=old2new_mapno)

    canonical_smi = reac_canonical_smi + ">>" + prod_canonical_smi
    canonical_smi_no_stereo = (
        reac_canonical_smi_no_stereo + ">>" + prod_canonical_smi_no_stereo
    )
    mol_no_stereo = (reac_mol_no_stereo, prod_mol_no_stereo)
    mol = (reac_mol, prod_mol)

    return canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol


def read_in_smiles(smiles):
    """
    Computes canonical SMILES strings and RDKit molecules from a list of SMILES strings.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES strings

    Returns
    -------
    canonical_smiles: List[str]
        List of canonical SMILES strings
    smiles_dict: dict
        Dictionary of canonical SMILES strings and their respective SMILES with stereoinformation and RDKit molecules.
    """

    canonical_smiles = []
    smiles_dict = {}
    for smi in smiles:
        check_if_reaction(smi, False)
        canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol, _ = preprocess(smi)
        smiles_dict[canonical_smi_no_stereo] = {
            "canonical_smi": [canonical_smi],
            "mol_diagram": mol_no_stereo,
            "mol_no_stereo": mol_no_stereo,
            "mol": [mol],
            "original_smi": smi,
        }
        canonical_smiles.append(canonical_smi_no_stereo)
    return canonical_smiles, smiles_dict


def read_in_smiles_unique(smiles):
    """
    Computes canonical SMILES strings and RDKit molecules from a list of SMILES strings and only keeps unique entries.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES strings

    Returns
    -------
    canonical_smiles: List[str]
        Sorted list of canonical SMILES strings
    smiles_dict: dict
        Dictionary of canonical SMILES strings and their respective SMILES with stereoinformation and RDKit molecules.
    skipped: List[str]
        List of skipped (because duplicate) SMILES strings.
    """

    canonical_smiles = []
    smiles_dict = {}
    skipped = []
    for smi in smiles:
        check_if_reaction(smi, False)
        canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol, _ = preprocess(smi)
        # New entry:
        if canonical_smi_no_stereo not in smiles_dict.keys():
            smiles_dict[canonical_smi_no_stereo] = {
                "canonical_smi": [canonical_smi],
                "mol_diagram": mol_no_stereo,
                "mol_no_stereo": mol_no_stereo,
                "mol": [mol],
                "original_smi": smi,
            }
            canonical_smiles.append(canonical_smi_no_stereo)
        # Known entry, but different stereochemistry:
        elif canonical_smi not in smiles_dict[canonical_smi_no_stereo]["canonical_smi"]:
            smiles_dict[canonical_smi_no_stereo]["canonical_smi"].append(canonical_smi)
            smiles_dict[canonical_smi_no_stereo]["mol"].append(mol)
        # Duplicate:
        else:
            skipped.append(smi)

    canonical_smiles = sorted(canonical_smiles)

    return canonical_smiles, smiles_dict, skipped


def read_in_reactions(rxn_smiles):
    """
    Computes canonical SMILES strings and RDKit molecules from a list of reaction SMILES strings.

    Parameters
    ----------
    rxn_smiles: List[str]
        List of reaction SMILES strings

    Returns
    -------
    canonical_smiles: List[str]
        List of canonical SMILES strings
    smiles_dict: dict
        Dictionary of canonical SMILES strings and their respective SMILES with stereoinformation and RDKit molecules.
    tags_core: dict
        Dictionary of atom map numbers of reaction center for each canonical SMILES.
    """

    canonical_smiles = []
    smiles_dict = {}
    tags_core = {}

    for rxn_smi in rxn_smiles:
        check_if_reaction(rxn_smi, True)
        canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol = preprocess_rxn(
            rxn_smi
        )
        change_ts_to_reac, change_ts_to_prod, ts, ts_seed = process_an_example(
            canonical_smi_no_stereo
        )
        tags = [atom.GetProp("molAtomMapNumber") for atom in ts_seed.GetAtoms()]

        smiles_dict[canonical_smi_no_stereo] = {
            "canonical_smi": [canonical_smi],
            "mol_diagram": ts,
            "mol_no_stereo": mol_no_stereo,
            "mol": [mol],
            "original_smi": rxn_smi,
        }
        tags_core[canonical_smi_no_stereo] = tags
        canonical_smiles.append(canonical_smi_no_stereo)

    return canonical_smiles, smiles_dict, tags_core


def read_in_reactions_unique(rxn_smiles):
    """
    Computes canonical SMILES strings and RDKit molecules from a list of reaction SMILES strings
    and only keeps unique entries.

    Parameters
    ----------
    rxn_smiles: List[str]
        List of reaction SMILES strings

    Returns
    -------
    seeds: List[str]
        List of SMILES seeds (of all minimal templates).
    rule_dict: dict
        A dictionary of all minimal templates of all seeds.
    canonical_smiles: List[str]
        List of canonical SMILES strings
    smiles_dict: dict
        Dictionary of canonical SMILES strings and their respective SMILES with stereoinformation and RDKit molecules.
    skipped: List[str]
        List of skipped (because duplicate) SMILES strings.
    tags_core: dict
        Dictionary of atom map numbers of reaction center for each canonical SMILES.
    """

    canonical_smiles = []
    smiles_dict = {}
    skipped = []
    tags_core = {}
    rule_dict = {}
    list_ts_seed = []
    seeds = []
    for rxn_smi in rxn_smiles:
        check_if_reaction(rxn_smi, True)
        canonical_smi_no_stereo, canonical_smi, mol_no_stereo, mol = preprocess_rxn(
            rxn_smi
        )
        change_ts_to_reac, change_ts_to_prod, ts, ts_seed = process_an_example(
            canonical_smi_no_stereo
        )
        tags = [atom.GetProp("molAtomMapNumber") for atom in ts_seed.GetAtoms()]
        for atom in ts_seed.GetAtoms():
            atom.SetNoImplicit(
                True
            )  # Set everything to "NoImplicit" in templates, we never want to add Hs without explicitly specifying
        ts_seed_smiles = moltosmiles_transition(
            ts_seed, {"reac": change_ts_to_reac, "prod": change_ts_to_prod}
        )

        # Create seed
        reac = Chem.MolFromSmiles(
            moltosmiles_transition(
                ts_seed, {"reac": change_ts_to_reac, "prod": change_ts_to_reac}
            ),
            sanitize=False,
        )
        prod = Chem.MolFromSmiles(
            moltosmiles_transition(
                ts_seed, {"reac": change_ts_to_prod, "prod": change_ts_to_prod}
            ),
            sanitize=False,
        )

        # New entry:
        if canonical_smi_no_stereo not in smiles_dict.keys():
            smiles_dict[canonical_smi_no_stereo] = {
                "canonical_smi": [canonical_smi],
                "mol_diagram": ts,
                "mol_no_stereo": mol_no_stereo,
                "mol": [mol],
                "original_smi": rxn_smi,
            }
            canonical_smiles.append(canonical_smi_no_stereo)
            tags_core[canonical_smi_no_stereo] = tags
        # Known entry, but different stereochemistry:
        elif canonical_smi not in smiles_dict[canonical_smi_no_stereo]["canonical_smi"]:
            smiles_dict[canonical_smi_no_stereo]["canonical_smi"].append(canonical_smi)
            smiles_dict[canonical_smi_no_stereo]["mol"].append(mol)
        # Duplicate:
        else:
            skipped.append(rxn_smi)
            break

        # Save seed if not new
        found_equal_seed = False
        for item in list_ts_seed:
            if (
                mols_are_equal(ts_seed, item[0])
                and mols_are_equal(reac, item[1])
                and mols_are_equal(prod, item[2])
            ):
                found_equal_seed = True
                rule_dict[item[3]]["smiles"].append(canonical_smi_no_stereo)
        if not found_equal_seed:
            seeds.append(ts_seed_smiles)
            rule_dict[ts_seed_smiles] = {
                "rule": ts_seed,
                "smiles": [canonical_smi_no_stereo],
                "change_dict": {"reac": change_ts_to_reac, "prod": change_ts_to_prod},
            }
            list_ts_seed.append((ts_seed, reac, prod, ts_seed_smiles))

    return seeds, rule_dict, canonical_smiles, smiles_dict, skipped, tags_core


def preprocess_seeds(seed_list, smiles, smiles_dict):
    """
    Preprocesses a list of seeds and matches them to a list of SMILES strings. If given an empty list of seeds,
    a seed with the maximum common substructure of all molecules is created.

    Parameters
    ----------
    seed_list: List[str]
        List of seeds.
    smiles: List[str]
        List of SMILES strings
    smiles_dict: dict
        Dictionary of canonical SMILES strings and their respective SMILES with stereoinformation and RDKit molecules.

    Returns
    -------
    seeds: List[str]
        List of seeds.
    rule_dict: dict
        A dictionary of all minimal templates of all seeds.
    num_smiles_seed: List[int]
        List of how many SMILES strings fit each seed.
    """

    if len(seed_list) == 0:
        mcs = find_common([smiles_dict[smi]["mol_no_stereo"] for smi in smiles])
        seeds = [Chem.MolToSmiles(mcs)]
    else:
        seeds = seed_list

    rule_dict = {}
    num_smiles_seed = []
    for seed in seeds:
        current_smiles = []
        current_rule = Chem.MolFromSmiles(seed, sanitize=False)
        for atom in current_rule.GetAtoms():
            atom.SetNoImplicit(
                True
            )  # Set everything to "NoImplicit" in templates, we never want to add Hs without explicitly specifying
        for smi in smiles:
            mol = smiles_dict[smi]["mol_no_stereo"]
            matches = mol.GetSubstructMatches(current_rule)
            for match in matches:
                if do_charges_fit(mol, current_rule, match):
                    current_smiles.append(smi)
                    break
        rule_dict[seed] = {
            "rule": current_rule,
            "smiles": current_smiles,
            "change_dict": {
                "reac": {"atom": {}, "bond": {}},
                "prod": {"atom": {}, "bond": {}},
            },
        }
        num_smiles_seed.append(len(current_smiles))

    return seeds, rule_dict, num_smiles_seed


def moltosmiles_transition(mol, change_dict):
    """
    Converts an RDKit mol with the changes specified in change_dict before converting to smiles.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    change_dict: dict
        A dictionary specifying changes to apply to bonds and atoms for reactants and products.

    Returns
    -------
    str
        Converted SMILES.
    """

    smiles = {}
    for mode in ["reac", "prod"]:
        mol2 = mol_transition(mol, change_dict[mode])
        smiles[mode] = Chem.MolToSmiles(mol2)

    if smiles["reac"] == smiles["prod"]:
        return smiles["reac"]
    else:
        return smiles["reac"] + ">>" + smiles["prod"]


def mol_transition(mol, change_dict):
    """
    Converts an RDKit mol with the changes specified in change_dict.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    change_dict: dict
        A dictionary specifying changes to apply to bonds and atoms.

    Returns
    -------
    mol2: rdkit.Chem.Mol
        Transformed RDKit molecule.
    """

    mol2 = deepcopy(mol)
    # Change bonds
    for idx in change_dict["bond"].keys():
        mol2.GetBondBetweenAtoms(idx[0], idx[1]).SetBondType(
            change_dict["bond"][idx][0]
        )
    # Delete bonds with zero bond order
    editable = Chem.EditableMol(mol2)
    for bond in reversed(mol2.GetBonds()):
        if bond.GetBondType() == Chem.rdchem.BondType.ZERO:
            editable.RemoveBond(
                bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            )
    mol2 = editable.GetMol()
    # Change charges
    for idx in change_dict["atom"].keys():
        mol2.GetAtomWithIdx(idx).SetFormalCharge(change_dict["atom"][idx][0])

    return mol2


def mols_are_equal(mol1, mol2):
    """
    Computes if two molecules are equal.

    Parameters
    ----------
    mol1: rdkit.Chem.Mol
        RDKit molecule.
    mol2: rdkit.Chem.Mol
        RDKit molecule.

    Returns
    -------
    bool
        Boolean whether molecules are equal.
    """

    if mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1):
        return True
    else:
        return False


def find_matching_atoms(train_mode, mol, rule, rule_small, tags_core):
    """
    Find the substructure matching of a rule or lowest matching rule to a molecule and confirm
    the match contains the reaction center.

    Parameters
    ----------
    train_mode: Literal["single_reactant", "transition_state"]
        Mode in which diagram was constructed.
    mol: rdkit.Chem.Mol
        RDKit molecule.
    rule: rdkit.Chem.Mol
        RDKit molecule of substructure/reaction rule.
    rule_small: rdkit.Chem.Mol
        RDKit molecule of substructure/reaction rule of lowest matching template
    tags_core: List[str]
        A list of the atom map numbers in the reaction center.

    Returns
    -------
    to_do_matches: list
        List of valid matches.
    """

    atoms_core = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.HasProp("molAtomMapNumber")
        and atom.GetProp("molAtomMapNumber") in tags_core
    ]
    matches = mol.GetSubstructMatches(rule)  # There might be multiple matches
    to_do_matches = []
    if rule_small:
        matches_small = mol.GetSubstructMatches(rule_small)
        for j in range(len(matches)):
            match_large = matches[j]
            if not match_includes_reaction_center(train_mode, match_large, atoms_core):
                continue
            for k in range(len(matches_small)):
                match = matches_small[k]
                if False in [x in match_large for x in match]:
                    continue
                if not match_includes_reaction_center(train_mode, match, atoms_core):
                    continue
                to_do_matches.append(match)
    else:
        for j in range(len(matches)):
            match = matches[j]
            if not match_includes_reaction_center(train_mode, match, atoms_core):
                continue
            to_do_matches.append(match)

    return to_do_matches


def match_includes_reaction_center(train_mode, match, atoms_core):
    """
    Determindes whether a substructure match includes the full reaction center.

    Parameters
    ----------
    train_mode: Literal["single_reactant", "transition_state"]
        Mode in which diagram was constructed.
    match: tuple
        Indices of substructure match.
    atoms_core: List[int]
        Atom indices belonging to the reaction center.

    Returns
    -------
    includes_rc: bool
        Boolean whether match includes the reaction center.
    """

    includes_rc = True
    if train_mode == "transition_state":
        if False in [core_atom in match for core_atom in atoms_core]:
            includes_rc = False
    return includes_rc


def make_fragment_indices(rule):
    """
    Fragments a molecule into multiple molecules where it is disconnected.
    For examples, turns 'CC.C' into 'CC' and 'C'.

    Parameters
    ----------
    rule: rdkit.Chem.Mol
        RDKit molecule.

    Returns
    -------
    rule_fragment_indices: tuple
        Tuple of tuple of atom indices in each fragment.
    rule_fragment_mols: tuple
        Tuple of RDKit molecule objects of each fragment.
    """

    rule_fragment_indices = Chem.GetMolFrags(rule)
    rule_fragment_mols = Chem.GetMolFrags(rule, asMols=True, sanitizeFrags=False)

    return rule_fragment_indices, rule_fragment_mols


def match_and_mol_transition(mol, change_dict, rule, direction="regular"):
    """
    Converts an RDKit mol with the changes specified in change_dict, with bond numbers corresponding to rule not mol
    thus an additional conversion step is needed to identify the correct bonds.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    change_dict: dict
        Changes to be applied to bonds and atoms, indexing corresponding to rule.
    rule: rdkit.Chem.Mol
        RDKit molecule of substructure/reaction rule.
    direction: Literal['regular', 'reversed'], default 'regular'
        Direction of change ('regular' or 'reversed').

    Returns
    mol2s: List[rdkit.Chem.Mol]
        List of transformed molecules.
    """

    matches = mol.GetSubstructMatches(rule)
    mol2s = []
    for match in matches:
        mol2 = deepcopy(mol)
        # Change bonds
        for idx in change_dict["bond"].keys():
            begin_atom_mol = match[idx[0]]
            end_atom_mol = match[idx[1]]
            bond = mol2.GetBondBetweenAtoms(begin_atom_mol, end_atom_mol)
            if direction == "regular":
                new_type = change_dict["bond"][idx][0]
            else:
                new_type = change_dict["bond"][idx][1]
            if bond is not None:
                bond.SetBondType(new_type)
            else:
                editable = Chem.EditableMol(mol2)
                editable.AddBond(begin_atom_mol, end_atom_mol, order=new_type)
                mol2 = editable.GetMol()

        # Delete bonds with zero bond order
        editable = Chem.EditableMol(mol2)
        for bond in reversed(mol2.GetBonds()):
            if bond.GetBondType() == Chem.rdchem.BondType.ZERO:
                editable.RemoveBond(
                    bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
                )
        mol2 = editable.GetMol()
        # Change charges
        for idx in change_dict["atom"].keys():
            if direction == "regular":
                mol2.GetAtomWithIdx(match[idx]).SetFormalCharge(
                    change_dict["atom"][idx][0]
                )
            else:
                mol2.GetAtomWithIdx(match[idx]).SetFormalCharge(
                    change_dict["atom"][idx][1]
                )
        mol2s.append(mol2)
    return mol2s


def mol_to_rxn_smiles(initial_smiles, initial_mols, d, predict_mode, verbose):
    """
    Creates a list of all possible reaction smiles from a given list of reactant(s).

    Parameters
    ----------
    initial_smiles: List[str]
        List of SMILES strings for reactant(s).
    initial_mols: List[rdkit.Chem.Mol]
        RDKit molecules of reactant(s) list.
    d: ehreact.diagram.diagram.Diagram
        Hasse Diagram.
    predict_mode: Literal["single_reactant", "multi_reactant", "transition_state"]
        Mode of prediction.
    verbose: bool
        Whether to print additional information.

    Returns
    -------
    current_smiles: List[str]
        List of reaction SMILES.
    belongs_to: List[int]
        List of indices which item in current_smiles belongs to which initial_smiles.
    combination: List[str]
        List of combination of reactants.
    """

    # Check for partners, create list of possibilities
    current_smiles = []
    belongs_to = []
    combination = []
    fragment_dict_reac = d.nodes[""].fragment_reac
    for i in range(len(initial_mols)):

        for key1 in fragment_dict_reac.keys():
            change_dict = d.nodes[""].change_dict[key1]
            rule_reac = d.nodes[key1].rule_reac
            rule = d.nodes[key1].rule

            # Check for partners, create list of possibilities
            if predict_mode == "single_reactant":
                missing_fragments = []
                for key2 in fragment_dict_reac[key1].keys():
                    rule_fragment = fragment_dict_reac[key1][key2]["rule_fragment"]
                    if not initial_mols[i].HasSubstructMatch(rule_fragment):
                        missing_fragments.append(
                            fragment_dict_reac[key1][key2]["molecules"]
                        )
                    else:
                        missing_fragments.append([initial_smiles[i]])
                combination_fragments = [
                    ".".join(tups)
                    for tups in list(itertools.product(*missing_fragments))
                ]

            # Simply use inputted molecules for multi_reactant mode
            else:
                combination_fragments = [initial_smiles[i]]

            combination_fragments = sorted(
                list(set(combination_fragments))
            )  # Only use unique entries

            # For each combination of reactants, create atom map numbers, update smiles, and create
            # the transition state and products according to change_dict and rule
            for combination_fragment in combination_fragments:
                if verbose:
                    print("Current fragment:", combination_fragment)
                mol = make_mol(combination_fragment)
                for atom in mol.GetAtoms():
                    atom.SetProp("molAtomMapNumber", str(atom.GetIdx() + 1))
                smiles_reac = Chem.MolToSmiles(mol)
                if verbose:
                    print("... Atom mapped reactants:", smiles_reac)
                ts_mols = match_and_mol_transition(
                    mol, change_dict["reac"], rule_reac, direction="reverse"
                )

                for ts_mol in ts_mols:
                    prod = match_and_mol_transition(
                        ts_mol, change_dict["prod"], rule, direction="regular"
                    )[0]
                    smiles_prod = Chem.MolToSmiles(prod)
                    if verbose:
                        print("... Possible product:", smiles_prod)
                    rxn_smiles = smiles_reac + ">>" + smiles_prod

                    # Check if valid rxn_smiles was produced (i.e. something changes):
                    if ".".join(sorted(smiles_reac.split("."))) == ".".join(
                        sorted(smiles_prod.split("."))
                    ):
                        continue
                    current_smiles.append(rxn_smiles)
                    belongs_to.append(i)
                    combination.append(
                        combination_fragment.replace(
                            initial_smiles[i] + ".", ""
                        ).replace("." + initial_smiles[i], "")
                    )

    return current_smiles, belongs_to, combination


def do_charges_fit(mol, current_rule, match):
    """
    Determines whether a current rule matches a molecule if formal charges are taken into account.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule object.
    current_rule: rdkit.Chem.Mol
        RDKit molecule object of rule
    match: tuple
         Indices of matching atoms in current_rule and mol.

    Returns
    -------
    do_charges_fit: bool
        Whether formal charges are the same.
    """

    do_charges_fit = True
    for idx in range(len(match)):
        if (
            current_rule.GetAtomWithIdx(idx).GetFormalCharge()
            != mol.GetAtomWithIdx(match[idx]).GetFormalCharge()
        ):
            do_charges_fit = False
    return do_charges_fit


def force_charge_fit(mol, current_rule, match):
    """
    Forces the formal charges of a rule to match the formal charges of a molecule.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule object.
    current_rule: rdkit.Chem.Mol
        RDKit molecule object of rule
    match: tuple
         Indices of matching atoms in current_rule and mol.

    Returns
    -------
    current_rule: rdkit.Chem.Mol
        RDKit molecule object of rule with updated formal charges.
    """

    for idx in range(len(match)):
        current_rule.GetAtomWithIdx(idx).SetFormalCharge(
            mol.GetAtomWithIdx(match[idx]).GetFormalCharge()
        )
    return current_rule
