from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

USE_IMAGE = True
# TODO use additional option to print tree to the screen (in case a user has no graphviz installed)


def draw_mol(mol, path, change_dict):
    """
    Function to transform and draw a smiles string using rdkit.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule.
    path: str
        Path to which a PNG image is saved.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    """

    try:
        AllChem.Compute2DCoords(mol)
        d = rdMolDraw2D.MolDraw2DCairo(200, 200)
        d.drawOptions().prepareMolsBeforeDrawing = False
        color_atoms = []
        color_bonds = []
        color_by_atoms = {}
        color_by_bonds = {}
        color = (0.8, 0.8, 0.8)
        color2 = (0.2, 0.8, 0.2)
        color3 = (0.8, 0.1, 0.20)
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")
                atom_idx = atom.GetIdx()
                color_atoms.append(atom_idx)
                color_by_atoms[atom_idx] = color

        for idx in change_dict["reac"]["atom"].keys():
            if (
                change_dict["prod"]["atom"][idx][0]
                > change_dict["reac"]["atom"][idx][0]
            ):
                color_by_atoms[idx] = color2
            elif (
                change_dict["prod"]["atom"][idx][0]
                < change_dict["reac"]["atom"][idx][0]
            ):
                color_by_atoms[idx] = color3
            else:
                color_by_atoms[idx] = color

        for idx in change_dict["reac"]["bond"].keys():
            bond = mol.GetBondBetweenAtoms(idx[0], idx[1])
            bond_idx = bond.GetIdx()
            color_bonds.append(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DATIVEL:
                color_by_bonds[bond_idx] = color2
                bond.SetBondType(change_dict["reac"]["bond"][idx][0])
            elif bond.GetBondType() == Chem.rdchem.BondType.DATIVER:
                color_by_bonds[bond_idx] = color3
                bond.SetBondType(change_dict["reac"]["bond"][idx][0])
            else:
                color_by_bonds[bond_idx] = color
            bond.SetBondType(change_dict["reac"]["bond"][idx][0])

        mol.UpdatePropertyCache()
        for atom in mol.GetAtoms():
            atom.GetExplicitValence()

        try:
            mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
        except:
            mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)

        if color_atoms != [] or color_bonds != []:
            d.DrawMolecule(
                mc,
                highlightAtoms=color_atoms,
                highlightBonds=color_bonds,
                highlightAtomColors=color_by_atoms,
                highlightBondColors=color_by_bonds,
            )
        else:
            d.DrawMolecule(mc)

        d.FinishDrawing()
        pic = d.GetDrawingText()
        with open(path, "wb") as f:
            f.write(pic)
    except Exception as e:
        print(e)
        print("Could not plot one of the molecules")


def write_beginning_dotfile(temp_dir_img):
    """
    Write the beginning of the dotfile (graphviz input).

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    """

    file = open(temp_dir_img + "/image_specs.txt", "w")

    file.write("digraph {\n")
    file.write("fontsize = 20\n")
    file.write("rank = sink\n")
    file.write("rankdir = RL\n")
    file.write('ranksep = "1.5"\n')
    file.write('nodesep = "0.2"\n')
    file.write('size = "100,100"')

    file.write("node[")
    file.write("penwidth = 5\n")
    file.write("color = blue\n")
    file.write('fillcolor = "#ffffff"\n')
    file.write("fontcolor = black\n")
    file.write("fontsize = 20\n")
    file.write("fontname = Helvetica\n")
    file.write("shape = square\n")
    file.write("style = filled\n]")

    file.write("edge [ arrowhead = open\n")
    file.write("color = black\n")
    file.write("fontcolor = black\n")
    file.write("fontname = Courier\n")
    file.write("fontsize = 12\n")
    file.write("style = solid\n")
    file.write("]\n")

    file.close()


def write_end_dotfile(temp_dir_img):
    """
    Write the end of the dotfile (graphviz input).

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    """

    file = open(temp_dir_img + "/image_specs.txt", "a")
    file.write("}")
    file.close()


def run_dot(temp_dir_img, save_plot):
    """
    Runs graphviz to plot the Hasse diagram.

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    save_plot: str
        File name to save PNG image to.
    """

    subprocess.check_call(
        "dot -o"
        + temp_dir_img
        + "/diagram.svg -Tsvg "
        + temp_dir_img
        + "/image_specs.txt",
        shell=True,
    )
    subprocess.check_call(
        "rsvg-convert " + temp_dir_img + "/diagram.svg >" + save_plot, shell=True
    )


def write_node(temp_dir_img, smiles, mol, is_leaf, file_numbers, change_dict):
    """
    Function to write a single node information to the image specification file

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    smiles: str
        SMILES string of the current node.
    mol: rdkit.Chem.Mol
        RDKit molecule to be drawn.
    is_leaf: bool
        Whether the current node is a leaf node.
    file_numbers: int
        Consecutive file identifier.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    """

    file = open(temp_dir_img + "/image_specs.txt", "a")
    file.write(double_quoted(smiles))
    file.write("[")
    file.write("label=" + double_quoted(""))
    if is_leaf:
        file.write("color=" + "black" + " ")
    else:
        file.write("color=" + "red" + " ")

    if USE_IMAGE and smiles != "":
        file.write(
            'image="' + temp_dir_img + "/" + file_numbers[smiles] + ".png" + '"\n'
        )
        draw_mol(
            mol=mol,
            path=temp_dir_img + "/" + file_numbers[smiles] + ".png",
            change_dict=change_dict,
        )
    file.write("shape=square ")
    file.write("]\n")


def write_edge(temp_dir_img, smiles1, smiles2):
    """
    Function to write a single edge information to the image specification file

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    smiles1: str
        SMILES string of the parent node.
    smiles2: str
        SMILES string of the child node.
    """

    file = open(temp_dir_img + "/image_specs.txt", "a")

    assert not ('"' in smiles1 or '"' in smiles2)

    file.write(double_quoted(smiles2))
    file.write(" -> ")
    file.write(double_quoted(smiles1))

    file.write("[dir=back ")
    file.write("shape=vee ")

    file.write("]\n")


def double_quoted(s):
    """
    Function to put double quotes around a string.

    Parameters
    ----------
    s: str
        String.

    Returns
    -------
    str
        String in double quotes.
    """

    return '"' + s + '"'


def remove_linear_chains(d):
    """
    Function to remove nodes that have exactly one parent and one child.

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
        A Hasse diagram.

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The modified Hasse diagram without linear chain nodes.
    """

    node_list = list(d.nodes)
    for n in node_list:
        key_node = n
        if key_node == "":
            continue
        key_parent = d.nodes[n].edges_to_parent[0].parent_node.key
        if len(d.nodes[n].edges_to_child) == 1:
            d.move_node(
                key_node=d.nodes[n].edges_to_child[0].child_node.key,
                key_parent_new=key_parent,
            )
            d.delete_edges_around_node(d.nodes[n].key)
    return d


def plot_hasse(d, temp_dir_img, save_plot, plot_only_branches):
    """
    Function to plot a Hasse diagram.

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
        A Hasse diagram.
    temp_dir_img: str
        Directory to save images to temporarily.
    save_plot: str
        File name to save PNG image to.
    plot_only_branches: bool
        Whether to remove linear chain nodes.
    """

    file_numbers = {}
    write_beginning_dotfile(temp_dir_img)
    smiles = ""
    seed = ""
    is_leaf = d.nodes[smiles].is_leaf
    mol = d.nodes[smiles].rule

    if plot_only_branches:
        d_compressed = remove_linear_chains(d)
        plot_iteration(
            temp_dir_img,
            d=d_compressed,
            smiles=smiles,
            mol=mol,
            is_leaf=is_leaf,
            file_numbers=file_numbers,
            seed=seed,
        )
    else:
        plot_iteration(
            temp_dir_img,
            d=d,
            smiles=smiles,
            mol=mol,
            is_leaf=is_leaf,
            file_numbers=file_numbers,
            seed=seed,
        )
    write_end_dotfile(temp_dir_img)
    run_dot(temp_dir_img, save_plot)


def plot_iteration(temp_dir_img, d, smiles, mol, is_leaf, file_numbers, seed):
    """
    Recursive function to plot the diagram. Per call, all current child nodes and corresponding edges are plotted.

    Parameters
    ----------
    temp_dir_img: str
        Directory to save images to temporarily.
    d: ehreact.diagram.diagram.Diagram
        A Hasse diagram.
    smiles: str
        SMILES string of the current node.
    mol: rdkit.Chem.Mol
        RDKit molecule to be drawn.
    is_leaf: bool
        Whether the current node is a leaf node.
    file_numbers: int
        Consecutive file identifier.
    seed: int
        Name of the minimal template of the current branch.
    """

    file_numbers[smiles] = str(len(file_numbers.keys()))
    if smiles != "" and d.nodes[""].change_dict != {}:
        change_dict = d.nodes[""].change_dict[seed]
    else:
        change_dict = {
            "reac": {"atom": {}, "bond": {}},
            "prod": {"atom": {}, "bond": {}},
        }

    write_node(
        temp_dir_img,
        smiles=smiles,
        mol=mol,
        is_leaf=is_leaf,
        file_numbers=file_numbers,
        change_dict=change_dict,
    )
    for edge in d.nodes[smiles].edges_to_child:
        child = edge.child_node.key
        write_edge(temp_dir_img, smiles2=child, smiles1=smiles)
        plot_iteration(
            temp_dir_img,
            d=d,
            smiles=child,
            mol=d.nodes[child].rule,
            is_leaf=d.nodes[child].is_leaf,
            file_numbers=file_numbers,
            seed=d.nodes[child].lowest_template,
        )
