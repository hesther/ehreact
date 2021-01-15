"""
arguments.py
A Python package for extracting and scoring reaction templates based on extended Hasse diagrams.

Handles the arguments of training and scoring Hasse diagrams.
"""

from tempfile import TemporaryDirectory
from tap import Tap  # pip install typed-argument-parser (https://github.com/swansonk14/typed-argument-parser)
from typing import List
from typing_extensions import Literal
import json
import os

class CommonArgs(Tap):
    """
    CommonArgs contains arguments that are used in both TrainArgs and PredictArgs.

    Attributes
    ----------
    no_qm: bool, default False
        Use QM information in building up the tree (precomputation of values) and during prediction in the scoring routine.
    verbose: bool, default False
        Whether to print additional information.
    quiet: bool, default False
        Whether to silence all output.
    stereochemistry: bool, default False
        Whether to use stereochemistry for scoring.
    compute_aam: bool, default False
        Whether to compute atom-mappings for reactions.
    """

    no_qm: bool = False
    verbose: bool = False
    quiet: bool = False
    stereochemistry: bool = False
    compute_aam: bool = False
    
    def process_args(self):
        if self.verbose is True and self.quiet is True:
            raise ValueError("Cannot use verbose and quiet output at the same time")
        if self.verbose is True:
            print("Outputting more information to screen")
        if self.stereochemistry:
            raise ValueError("Currently not implemented, please have some patience. Exit program.")


class TrainArgs(CommonArgs):
    """
    TrainArgs includes CommonArgs along with additional arguments used for generating a template tree from a list of reaction smiles.

    Attributes
    ----------
    data_path: str
        Path to data CSV file.
    save_path: str, default  None
        File to which diagram is saved.
    plot:  bool, default False
        Boolean whether the diagram is plotted.
    train_mode: Literal["single_reactant","transition_state"], default "transition_state" 
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.
    seed: List[str], default [] 
        List of SMILES seeds for the reactant algorithm, usually a single seed is given.
    negative_data: str, default None 
        Path to CSV file with negative data.
    no_props: bool, default False 
        Do not compute any properties, just output the diagram.
    plot_only_branches: bool, default False 
        Plot only substructures that branch off.
    """

    data_path: str 
    save_path: str = None 
    plot:  bool = False
    train_mode: Literal["single_reactant","transition_state"] = "transition_state" 
    seed: List[str] = []
    negative_data: str = None 
    no_props: bool = False 
    plot_only_branches: bool = False

    def process_args(self):
        super(TrainArgs, self).process_args()

        global temp_dir  # Prevents the temporary directory from being deleted upon function return

        if self.plot: #Create temporary directory for plotting routine
            temp_dir = TemporaryDirectory()
            self.temp_dir_img = temp_dir.name
                        
        if self.seed != [] and self.train_mode != "single_reactant":
            raise ValueError("A seed can only be specified when working in reactant mode"
                             "Please specify --train_mode single_reactant or omit the --seed.")

        if self.seed == [] and self.train_mode == 'single_reactant':
            raise ValueError("A seed must be specified when working in reactant mode."
                             "Specify '' (empty string) if you want to autogenerate a MCS seed"
                             "Specify --seed <seed>, e.g. --seed 'C([H])O[H]' or multiple (exclusive) seeds as  --seed 'C([H])O[H]' 'C=O'")

        if self.compute_aam and self.train_mode is not 'transition_state':
            raise ValueError("Computing atom maps is only possible for reactions, please specify --train_mode transition_state.")

class PredictArgs(CommonArgs):
    """
    PredictArgs includes CommonArgs along with additional arguments used for predicting reaction scores from a template tree and a list of reaction smiles.

    Attributes
    ----------
    load_path: str, default None 
        File from which diagram is loaded.
    test_path: str, default None
        Path to CSV file containing testing data for which predictions will be made.
    preds_path: str, default None 
        Path to CSV file where predictions will be saved.
    predict_mode: Literal["single_reactant","multi_reactant","transition_state"], default "transition_state" 
        Predict mode, either transition states (input reaction smiles, with train_mode='transition_state') or multiple reactants (input list of smiles, with train_mode='transition_state') or single reactants (input smiles, compatible with both train_modes).
    hyper_params: str, default None 
        File from which hyperparameters are loaded.
    save_rawnumber: str, default None 
        File to which raw scores can be saved, of interest for fast hyperparameter tuning.
    """

    load_path: str = None
    test_path: str = None 
    preds_path: str = None
    predict_mode: Literal["single_reactant","multi_reactant","transition_state"] = "transition_state" 
    hyper_params: str = None 
    save_rawnumber: str = None 

    def load_params(self):
        with open(self.hyper_params) as json_file:
            params = json.load(json_file)
        return params

    def process_args(self):
        super(PredictArgs, self).process_args()

        if self.hyper_params:
            self.params=self.load_params()
        else:
            self.params={}
            self.params["hyper_template"]={"use_template":"lowest"}
            self.params["hyper_location"]={"n":2}
            self.params["hyper_diversity"]={"mode":"mean","shift_amount":0.5,"cap":1,"cap_mean":0.8,"cap_std":0.05,"cap_via_overall":True}
            self.params["hyper_overall"]={"mode":"sigmoid","center":2.5,"shift":1}
            self.params["hyper_similarity"]={"mode":"mae","factors":[["partial_charge",3]],"use_prod":True}
            self.params["hyper_negative"]= {"use": True,"n_sim": 2}
            self.params["hyper_shape"]={"use": False,"charges":"qm_pred"}
            self.params["hyper_score"]={"mode":"mult","factors":[1,1,1,1,1,1,1]}
        if not self.quiet:
            print(self.params)

        if self.compute_aam and self.predict_mode is not 'transition_state':
            raise ValueError("Computing atom maps is only possible for reactions, please specify --predict_mode transition_state.")

