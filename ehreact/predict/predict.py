"""
Entry point for predicting on a calculated Hasse diagram.
"""
import csv
import pickle
from ehreact.predict import make_prediction


def predict(args):
    """
    Scores queries on a precomputed Hasse diagram based on the inputted arguments

    Parameters
    ----------
    args: Namespace
         Namespace of arguments.

    """

    if not args.quiet:
        print(args)

    # Read in test data:
    with open(args.test_path) as f:
        smiles = f.read().splitlines()

    # Load diagram:
    with open(args.load_path, "rb") as handle:
        d = pickle.load(handle)

    scores, combination, current_smiles, belongs_to, _ = make_prediction(
        smiles=smiles,
        d=d,
        params=args.params,
        verbose=args.verbose,
        quiet=args.quiet,
        compute_aam=args.compute_aam,
        predict_mode=args.predict_mode,
        stereochemistry=args.stereochemistry,
    )

    if args.preds_path:
        with open(args.preds_path, "w") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "Input SMILES",
                    "Score",
                    "Additional reaction partner",
                    "Scored reaction",
                ]
            )
            for i in range(len(current_smiles)):
                writer.writerow(
                    [
                        smiles[belongs_to[i]],
                        scores[i],
                        combination[i],
                        current_smiles[i],
                    ]
                )
        if not args.quiet:
            print("Saved predictions to", args.preds_path)
