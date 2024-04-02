#!/usr/bin/python3.8

from __future__ import annotations

import contextlib
import functools
import glob
import itertools
import multiprocessing
import os
import pathlib
import subprocess
import sys

print = functools.partial(print, flush=True)


@functools.singledispatch
def align_reads() -> None:
    raise NotImplementedError


@align_reads.register
def _(read_file: str, genome: str) -> None:
    ogd = os.getcwd()
    directory, _, file = read_file.rpartition('/')
    os.chdir(directory)

    try:
        print(f'Aligning {read_file=} in {os.getcwd()}')
        prefix = file.partition('.')[0]

        subprocess.run(
            [
                'STAR',
                '--runMode', 'alignReads',
                '--runThreadN', '8',
                '--outFilterMultimapNmax', '20',
                '--outFileNamePrefix', prefix,
                '--genomeDir', genome,
                '--readFilesIn', file,
                '--readFilesCommand', 'zcat' if read_file.endswith('gz') else 'cat',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'MD'
            ],
            stdin=None, stdout=sys.stdout, stderr=sys.stderr,
            check=True
        )

        print(f'Finished aligning {read_file=}')
    finally:
        with contextlib.suppress(OSError):
            os.chdir(ogd)


@align_reads.register
def _(directories: list, genome: str, processes: int = 8) -> None:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            align_reads,
            zip(
                itertools.chain.from_iterable(
                    glob.glob(f'{directory}/**/*.fq.gz', recursive=True) for directory in directories
                ),
                itertools.repeat(genome)
            )
        )
        

def main():
    import argparse

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('directories', nargs='*', default='.', type=str)
    parser.add_argument('-g', '--genome', type=str)
    parser.add_argument('-p', '--processes', default=8, type=int)

    args = parser.parse_args()
    print(vars(args))

    print('---BEGIN---')
    align_reads(args.directories, str(pathlib.Path(args.genome).absolute()), args.processes)
    print('---END---')


if __name__ == '__main__':
    main()
