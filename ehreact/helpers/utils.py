import itertools


def findsubsets(S, m):
    """
    Find all subsets of S containing m elements.

    Parameters
    ----------
    S: list
        List of objects.
    m: int
        Number of elements per subset.

    Returns
    -------
    subset_list: list
        List of subsets.
    """

    subset = list(itertools.combinations(S, m))
    subset_list = [list(x) for x in subset]
    return subset_list
