"""
arguments.py
Handles the arguments of training and scoring Hasse diagrams.
"""

from tempfile import TemporaryDirectory
from tap import Tap
from typing import List
from typing_extensions import Literal
import json


class CommonArgs(Tap):
    """
    CommonArgs contains arguments that are used in both TrainArgs and PredictArgs.

    Attributes
    ----------
    verbose: bool, default False
        Whether to print additional information.
    quiet: bool, default False
        Whether to silence all output.
    compute_aam: bool, default False
        Whether to compute atom-mappings for reactions.
    """

    verbose: bool = False
    quiet: bool = False
    compute_aam: bool = False

    def process_args(self):
        if self.verbose is True and self.quiet is True:
            raise ValueError("Cannot use verbose and quiet output at the same time")
        if self.verbose is True:
            print("Outputting more information to screen")


class TrainArgs(CommonArgs):
    """
    TrainArgs includes CommonArgs along with additional arguments used for generating a
    template tree from a list of reaction smiles.

    Attributes
    ----------
    data_path: str
        Path to data CSV file.
    save_path: str, default  None
        File to which diagram is saved.
    save_plot:  str, default None
        File to which save image of diagram.
    train_mode: Literal["single_reactant", "transition_state"], default "transition_state"
        Train mode, either transition states extracted from reaction smiles or single reactants extracted from smiles.
    seed: List[str], default []
        List of SMILES seeds for the reactant algorithm, usually a single seed is given.
    no_props: bool, default False
        Do not compute any properties, just output the diagram.
    plot_only_branches: bool, default False
        Plot only substructures that branch off.
    """

    data_path: str
    save_path: str = None
    save_plot: str = None
    train_mode: Literal["single_reactant", "transition_state"] = "transition_state"
    seed: List[str] = []
    no_props: bool = False
    plot_only_branches: bool = False

    def process_args(self):
        super(TrainArgs, self).process_args()

        global temp_dir  # Prevents the temporary directory from being deleted upon function return

        if self.save_plot:  # Create temporary directory for plotting routine
            temp_dir = TemporaryDirectory()
            self.temp_dir_img = temp_dir.name
        else:
            self.temp_dir_img = None

        if self.seed != [] and self.train_mode != "single_reactant":
            raise ValueError(
                "A seed can only be specified when working in reactant mode "
                "Please specify --train_mode single_reactant or omit the --seed."
            )

        if self.compute_aam and self.train_mode != "transition_state":
            raise ValueError(
                "Computing atom maps is only possible for reactions, please specify --train_mode transition_state."
            )


class PredictArgs(CommonArgs):
    """
    PredictArgs includes CommonArgs along with additional arguments used for predicting reaction scores
    from a template tree and a list of reaction smiles.

    Attributes
    ----------
    load_path: str
        File from which diagram is loaded.
    test_path: str
        Path to CSV file containing testing data for which predictions will be made.
    preds_path: str, default None
        Path to CSV file where predictions will be saved.
    predict_mode: Literal["single_reactant", "multi_reactant", "transition_state"], default "transition_state"
        Predict mode, either transition states (input reaction smiles, with train_mode='transition_state')
        or multiple reactants (input list of smiles, with train_mode='transition_state') or single reactants
        (input smiles, compatible with both train_modes).
    hyper_params: str, default None
        File from which hyperparameters are loaded.
    stereochemistry: bool, default False
        Whether to use stereochemistry for scoring.
    """

    load_path: str
    test_path: str
    preds_path: str = None
    predict_mode: Literal[
        "single_reactant", "multi_reactant", "transition_state"
    ] = "transition_state"
    hyper_params: str = None
    stereochemistry: bool = False

    def process_args(self):
        super(PredictArgs, self).process_args()

        if self.hyper_params:
            with open(self.hyper_params) as json_file:
                self.params = json.load(json_file)
        else:
            self.params = None

        if self.compute_aam and self.predict_mode != "transition_state":
            raise ValueError(
                "Computing atom maps is only possible for reactions, please specify --predict_mode transition_state."
            )
