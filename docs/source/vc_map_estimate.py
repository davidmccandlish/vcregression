import argparse
parser = argparse.ArgumentParser()
parser.add_argument("a", help="alphabet size", type=int)
parser.add_argument("l", help="sequence length", type=int)
parser.add_argument("-name", help="name of output folder")
parser.add_argument("-data", help="path to input data",  type=str, required=True)
parser.add_argument("-lambdas", help="path to lambdas",
                    type=str, required=True)

