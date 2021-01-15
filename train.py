"""Generates template trees on a dataset."""

from ehreact.arguments import TrainArgs
from ehreact.train import train

if __name__ == '__main__':
    args = TrainArgs().parse_args()
    train(args)
