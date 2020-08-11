import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-samples", help="tuple of paths to HMC samples",
                    dest="samplepaths", nargs="+", type=str, required=True)
parser.add_argument("-name", help="project name")


