#!/usr/bin/env python3

import argparse
from typing import Tuple

from charge.repository import Repository


def get_args() -> Tuple[str, str]:
    parser = argparse.ArgumentParser(
            description='Create a repository from a directory of LGF files.\n'
                        'This program reads molecule definitions from LGF'
                        ' files in a given directory, and processes them into'
                        ' a repository file for use with charge_assign.')
    parser.add_argument('-i', '--input', help='Path of the directory with input'
                        ' files')
    parser.add_argument('-r', '--repository', help='Path to the repository file'
                        ' to create')
    parser.add_argument('-m', '--max-shell', type=int, default=3,
                        help='Maximum shell size to use')
    args = parser.parse_args()
    return args.input, args.repository, args.max_shell


def main() -> None:
    source_dir, repo_file, max_shell = get_args()

    repo = Repository.create_from(source_dir, max_shell=max_shell)
    repo.write(repo_file)


if __name__ == '__main__':
    main()
