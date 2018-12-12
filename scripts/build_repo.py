#!/usr/bin/env python3

import argparse
from typing import Tuple

from charge.repository import Repository


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            description='Create a repository from a directory of LGF files.\n'
                        'This program reads molecule definitions from LGF'
                        ' files in a given directory, and processes them into'
                        ' a repository file for use with charge_assign.')
    parser.add_argument('-i', '--input', help='Path of the directory with input'
                        ' files', required=True)
    parser.add_argument('-r', '--repository', help='Path to the repository file'
                        ' to create', required=True)
    parser.add_argument('-m', '--max-shell', type=int, default=3,
                        help='Maximum shell size to use')
    parser.add_argument('-t', '--traceable', help='Generate a traceable'
                        ' repository', action='store_true')
    args = parser.parse_args()
    return args


def main() -> None:
    args = get_args()

    repo = Repository.create_from(args.input, max_shell=args.max_shell,
                                  traceable=args.traceable)
    repo.write(args.repository)


if __name__ == '__main__':
    main()
