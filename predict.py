"""Generates template trees on a dataset."""

from ehreact.arguments import PredictArgs
from ehreact.predict import predict

if __name__ == '__main__':
    args = PredictArgs().parse_args()
    predict(args)
