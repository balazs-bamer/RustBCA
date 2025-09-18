import os
import sys
from rustbcaPly import *


def main() -> int:
    scale = 0.0001
    if len(sys.argv) < 2:
        raise Exception(f'usage: {sys.argv[0]} filename [abstract scale = {scale}]')
    nameIn = sys.argv[1]
    if len(sys.argv) > 2:
        scale = float(sys.argv[2])
    nameBase = ''
    was = ''
    strOptions = '[options]'
    with open(nameIn, encoding='utf-8') as file:
        for line in file:
            tokens = line.split()
            if tokens[0] == strOptions:
                was = strOptions
            elif tokens[0] == 'name' and was == strOptions:
                nameBase = tokens[2].replace('"', '')
                break;
    do_trajectory_plot_3d(nameBase, input_file=nameIn, radius=0.1, scaleAbstr=scale)


if __name__ == '__main__':
    sys.exit(main())

