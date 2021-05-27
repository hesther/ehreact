"""
train.py
Entry point for training Hasse diagrams.
"""
from ehreact.train import calculate_diagram


def train(args):
    """
    Computes a Hasse diagram based on the inputted arguments

    Parameters
    ----------
    args: Namespace
         Namespace of arguments.
    """

    if not args.quiet:
        print(args)

    # Read in positive data:
    with open(args.data_path) as f:
        smiles = f.read().splitlines()

    _ = calculate_diagram(
        smiles=smiles,
        verbose=args.verbose,
        quiet=args.quiet,
        compute_aam=args.compute_aam,
        save_path=args.save_path,
        save_plot=args.save_plot,
        train_mode=args.train_mode,
        seed=args.seed,
        no_props=args.no_props,
        plot_only_branches=args.plot_only_branches,
        temp_dir_img=args.temp_dir_img,
    )
