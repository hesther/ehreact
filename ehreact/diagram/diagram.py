class Edge:
    """
    An Edge connecting the parent to the child node.

    Attributes
    ----------
    parent_node: ehreact.diagram.diagram.Node
        Parent node of the edge (more general template).
    child_node: ehreact.diagram.diagram.Node
        Child node of the edge (more specific template).
    """

    def __init__(self, parent_node, child_node):
        """
        Initializes and edge.

        Parameters
        ----------
        parent_node: ehreact.diagram.diagram.Node
            Parent node of the edge (more general template).
        child_node: ehreact.diagram.diagram.Node
            Child node of the edge (more specific template).
        """

        self.parent_node = parent_node
        self.child_node = child_node
        self.parent_node.edges_to_child.append(self)
        self.child_node.edges_to_parent.append(self)

    def __str__(self):
        """
        Sets the informal string representation of the Edge.

        Returns
        -------
        str
            String of the parent and child node key.
        """

        return str(self.parent_node) + "---" + str(self.child_node)

    def __repr__(self):
        """
        Sets the formal string representation of the Edge.

        Returns
        -------
        str
            String of the parent and child node key.
        """

        return str(self.parent_node) + "---" + str(self.child_node)


class Node:
    """
    A Node holding information to which parent and child(s) it relates to, which rule it was derived from
    and additional data.

    Attributes
    ----------
    edges_to_child: List[ehreact.diagram.diagram.Edge]
        List of edges to child nodes.
    edges_to_parent: listList[ehreact.diagram.diagram.Edge]
        List of edges to parent nodes.
    key: str
        Name of node.
    is_leaf: bool
        Whether the node is a leaf node (has no children).
    rule: rdkit.Chem.Mol
        RDKit molecule object of template (molecule or pseudomolecule).
    rule_reac: rdkit.Chem.Mol
        RDKit molecule object of reactant side of template (for transition state mode).
    rule_prod: rdkit.Chem.Mol
        RDKit molecule object of product side of template (for transition state mode).
    all_leafs: List[str]
        List of all SMILES strings of leaf nodes.
    min_dist_leafs: int
        Minimum distance to nearest leaf node.
    fp_reac: List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
        List of fingerprints of reactants without stereoinformation.
    fp_reac_stereo: List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
        List of fingerprints of reactants with stereoinformation.
    fp_prod: List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
        List of fingerprints of products without stereoinformation.
    fp_prod_stereo: List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
        List of fingerprints of products with stereoinformation.
    diversity_reac: float
        Mean pair similarity of all reactants of leaf nodes in the branch without stereoinformation.
    diversity_prod: float
        Mean pair similarity of all products of leaf nodes in the branch without stereoinformation.
    diversity_reac_stereo: float
        Mean pair similarity of all reactants of leaf nodes in the branch with stereoinformation.
    diversity_prod_stereo: float
        Mean pair similarity of all products of leaf nodes in the branch with stereoinformation.
    tags_core: List[str]
        A list of the atom map numbers in the reaction center.
    fragment_reac: dict
        A dictionary of fragments for a reaction.
    change_dict: dict
        A dictionary of all changes upon going from reactants to products.
    lowest_template: str
        Name of the lowest template of the current branch.
    """

    def __init__(self, key, rule, is_leaf):
        """
        Initializes a node.

        Parameters
        ----------
        key: str
            Unique name of the node, usually SMILES string.
        rule: rdkit.Chem.Mol
            RDKit molecule of the reaction rule.
        is_leaf: bool
            Boolean whether the node is a leaf node (full molecule or reaction) or a reaction rule
            (substructure or reaction template).
        """

        self.edges_to_child = []
        self.edges_to_parent = []
        self.key = key
        self.is_leaf = is_leaf
        self.rule = rule
        self.rule_reac = None
        self.rule_prod = None
        self.all_leafs = list()
        self.min_dist_leaf = 9999
        self.fp_reac = []
        self.fp_reac_stereo = []
        self.fp_prod = []
        self.fp_prod_stereo = []
        self.diversity_reac = None
        self.diversity_prod = None
        self.diversity_reac_stereo = None
        self.diversity_prod_stereo = None
        self.tags_core = None
        self.fragment_reac = None
        self.change_dict = {}
        self.lowest_template = None

    def __str__(self):
        """
        Sets the informal string representation of the Node.

        Returns
        -------
        str
            String of the node key.
        """

        return "root" if self.key == "" else self.key

    def __repr__(self):
        """
        Sets the formal string representation of the Node.

        Returns
        -------
        str
            String of the node key.
        """

        return "root" if self.key == "" else self.key

    @classmethod
    def create_root_node(cls):
        """
        Function to create an empty root node.

        Returns
        -------
        ehreact.diagram.diagram.Node
            An empty root node
        """

        n = Node(key="", rule=None, is_leaf=False)
        return n

    def create_new_node(key, rule, is_leaf):
        """
        Function to create a new node.

        Parameters
        ----------
        key: str
            Key of the newly created node.
        rule: rdkit.Chem.Mol
            RDKit molecule of the substructure/reaction rule of the node
        is_leaf: bool
            Boolean whether new node is a leaf node.

        Returns
        -------
        ehreact.diagram.diagram.Node
            A new node
        """

        return Node(key=key, rule=rule, is_leaf=is_leaf)


class Diagram:
    """
    Custom-made diagram class for extended Hasse diagrams.

    Attributes
    ----------
    nodes: dict
        Dictionary of all nodes.
    root: ehreact.diagram.diagram.Node
        Root node.
    mode: str
        Mode in which diagram was calculated.
    """

    def __init__(self):
        """
        Initializes a diagram.
        """

        self.nodes = dict()
        self.root = Node.create_root_node()
        self.nodes[self.root.key] = self.root
        self.mode = None

    def add_node(self, key_node, rule, is_leaf, key_parent_node, quiet):
        """
        Function to add a new node to the diagram.

        Parameters
        ----------
        key_node: str
            Key of the newly created node.
        rule: rdkit.Chem.Mol
            RDKit molecule of the substructure/reaction rule of the node
        is_leaf: bool
            Boolean whether new node is a leaf node.
        key_parent_node: str
            Key of the parent node of the newly created node.
        """

        n = Node.create_new_node(key=key_node, rule=rule, is_leaf=is_leaf)
        self.nodes[n.key] = n
        Edge(self.nodes[key_parent_node], n)
        if not quiet:
            print("...", self.nodes[key_parent_node], "to", n)

    def insert_node(self, key_node, rule, is_leaf, key_parent_node, key_child_node):
        """
        Function to insert a node between an existing parent-child pair.

        Parameters
        ----------
        key_node: str
            Key of the newly created node.
        rule: rdkit.Chem.Mol
            RDKit molecule of the substructure/reaction rule of the node
        is_leaf: bool
            Boolean whether new node is a leaf node.
        key_parent_node: str
            Key of the parent node of the newly created node.
        key_child_node: str
            Key of the child node of the newly created node.
        """

        n = Node.create_new_node(key=key_node, rule=rule, is_leaf=is_leaf)
        self.nodes[n.key] = n

        # Delete old edge
        del_edge = self.nodes[key_child_node].edges_to_parent[0]
        del_edge.child_node.edges_to_parent.remove(del_edge)
        del_edge.child_node = None
        del_edge.parent_node.edges_to_child.remove(del_edge)
        del_edge.parent_node = None

        # Form two new edges
        Edge(self.nodes[key_parent_node], n)
        Edge(n, self.nodes[key_child_node])

    def move_node(self, key_node, key_parent_new):
        """
        Function to attach a child node to a new parent node.

        Parameters
        ----------
        key_node: str
            Key of the node that is attached to a new parent.
        key_parent_new: str
            Key of the new parent.
        """

        # Delete old edge
        del_edge = self.nodes[key_node].edges_to_parent[0]
        del_edge.child_node.edges_to_parent.remove(del_edge)
        del_edge.child_node = None
        del_edge.parent_node.edges_to_child.remove(del_edge)
        del_edge.parent_node = None

        # Form new edge
        Edge(self.nodes[key_parent_new], self.nodes[key_node])

    def delete_edges_around_node(self, key_node):
        """
        Function to remove a node from the diagram, with removing all edges.
        The node will still be saved in the diagram, but without any edges connected to it.
        This is only meaningful for plotting, not a useful function elsewhere!

        Parameters
        ----------
        key_node: str
            Key of the node that is to be removed.
        """

        # Delete edge to parent
        del_edges = (
            self.nodes[key_node].edges_to_parent + self.nodes[key_node].edges_to_child
        )
        for del_edge in del_edges:
            del_edge.child_node.edges_to_parent.remove(del_edge)
            del_edge.child_node = None
            del_edge.parent_node.edges_to_child.remove(del_edge)
            del_edge.parent_node = None


def sanity_check(d):
    """
    Perform a sanity check on the diagram.

    Parameters
    ----------
    d: ehreact.diagram.diagram.Diagram
        Hasse Diagram.
    """

    # Diagram not empty:
    if len(d.nodes) <= 1:
        raise ValueError("Empty diagram")

    # Key is never used twice:
    if len(d.nodes.keys()) > len(set(d.nodes.keys())):
        raise ValueError("Key occured multiple times")

    # Every node (except root) has exactly one parent:
    for n in d.nodes:
        if d.nodes[n].key != "" and len(d.nodes[n].edges_to_parent) != 1:
            raise ValueError("Found a node with number of parent != 1")
