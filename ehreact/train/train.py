"""
Entry point for training Hasse diagrams.
"""

from argparse import Namespace
from ehreact.train import calculate_diagram

def train(args):
    """
    Computes a Hasse diagram based on the inputted arguments

    Parameters
    ----------
    args: Namespace
         Namespace of arguments.

    """

    #Preload QM machine learning model
    if not args.no_qm:
        from ehreact.preprocess.chemprop_qm_prediction import load_chemprop_model
        chemprop_args,train_args,scaler,features_scaler,model=load_chemprop_model()
        chemprop_args=Namespace(args=chemprop_args,train_args=train_args,scaler=scaler,features_scaler=features_scaler,model=model)
    else:
        chemprop_args=Namespace()

    #Read in positive data:
    with open(args.data_path) as f:
        smiles = f.read().splitlines()

    #Read in negative data:
    if args.negative_data != None:
        with open(args.negative_data) as f:
            negative_smiles = f.read().splitlines()
    else:
        negative_smiles=[]
        
    calculate_diagram(smiles=smiles,
                      negative_smiles=negative_smiles,
                      chemprop_args=chemprop_args,
                      no_qm=args.no_qm,
                      verbose=args.verbose,
                      quiet=args.quiet,
                      stereochemistry=args.stereochemistry,
                      compute_aam=args.compute_aam,
                      save_path=args.save_path,
                      plot=args.plot,
                      train_mode=args.train_mode,
                      seed=args.seed,
                      no_props=args.no_props,
                      plot_only_branches=args.plot_only_branches
    )
