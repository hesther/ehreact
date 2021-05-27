import numpy as np


def highest_sim(sims_reac, sims_prod, use_prod):
    """
    Function to find highest similarity of a query to a set of leaf nodes.

    Parameters
    ----------
    sims_reac: List[float]
        List of similarities of the reactants.
    sims_prod: List[float]
        List of similarities of the products.
    use_prod: bool
        Whether to use the product similarities in addition to the reactants.

    Returns
    -------
    best_score: float
        Highest similarity score.
    """

    best_score = 0
    for i in range(len(sims_reac)):
        if use_prod:
            best_score = max(best_score, np.mean([sims_reac[i], sims_prod[i]]))
        else:
            best_score = max(best_score, sims_reac[i])
    return best_score


def default_score_generator(rawnumber, params, verbose):
    """
    Score generator.

    Parameters
    ----------
    rawnumber: dict
        Dictionary of raw scores.
    params: dict
        Dictionary of hyperparameters.
    verbose: bool
        Whether to print additional information.

    Returns
    -------
    best_score: float
        Overall score.
    """

    if rawnumber["train_mode"] == "transition_state":
        hyper_similarity_use_prod = params["use_prod"]
    else:
        hyper_similarity_use_prod = False  # No product fps for single_reactant mode

    sc_sim = highest_sim(
        rawnumber["chemical_scores"],
        rawnumber["chemical_scores_prod"],
        hyper_similarity_use_prod,
    )
    if hyper_similarity_use_prod:
        sc_sim_ov = np.mean(
            [
                np.mean(rawnumber["chemical_scores_overall"]),
                np.mean(rawnumber["chemical_scores_overall_prod"]),
            ]
        )
        sc_div_ov = min(
            np.mean(
                [rawnumber["overall_div_mean_reac"], rawnumber["overall_div_mean_prod"]]
            ),
            params["cap_sp"],
        )
    else:
        sc_sim_ov = np.mean(rawnumber["chemical_scores_overall"])
        sc_div_ov = min(rawnumber["overall_div_mean_reac"], params["cap_sp"])
    sc_dist = min(rawnumber["min_dist_leaf"], params["cap_sl"]) - 1
    best_score = (
        params["coeff_ss"]
        * sc_sim
        * (
            1
            + params["coeff_sm"] * sc_sim_ov
            + params["coeff_sp"] * sc_div_ov
            + params["coeff_sl"] * sc_dist
        )
    )
    best_score = min(best_score, 1)
    best_score = max(best_score, 0)

    if verbose:
        print("S_S:", sc_sim, "S_P", sc_div_ov, "S_M", sc_sim_ov, "S_L", sc_dist)

    return best_score
