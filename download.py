#!/usr/bin/python3.8

from __future__ import annotations

import bs4
import concurrent.futures
import contextlib
import ftplib
import functools
import itertools
import pathlib
import re
import requests
import string


NCBI_ACCESSION = re.compile('^\S{3}\d{7}$')
NCBI_GEO_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
NCBI_SRA_URL = 'https://www.ncbi.nlm.nih.gov/sra'
ENA_URL = 'ftp.sra.ebi.ac.uk'

parse = functools.partial(bs4.BeautifulSoup, features='html.parser')
print = functools.partial(print, flush=True)


def search_accession(s: str) -> bool:
    return s and NCBI_ACCESSION.match(s)


def search_for(target: str):
    target = target.lower().strip()
    def search(s: str) -> bool:
        return s and s.lower().strip() == target

    return search


def download_samples(accessions: list[str]) -> None:
    def download(accession: str) -> None:
        read = pathlib.Path('Reads')
        read.mkdir(parents=True, exist_ok=True)

        with ftplib.FTP(ENA_URL) as ftp:
            ftp.login()
            
            print(f'Beginning download of {accession=}')
            ftp.cwd(f'/vol1/fastq/{accession[:6]}/00{accession[-1]}/{accession}')

            for file in ftp.nlst():
                if (read / file).exists():
                    print(f'{file=} already exists... moving on')
                    continue

                print(f'Beginning download of {file=}')
                tmp = read / f'tmp_{file}'
                try:
                    with tmp.open(mode='wb') as fp:
                        ftp.retrbinary(f'RETR {file}', fp.write)

                    tmp.rename(read / file)
                finally:
                    tmp.unlink(missing_ok=True)

                print(f'Finished download of {file=}')
            print(f'Finished download of {accession=}')

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for accession in accessions:
            pool.submit(download, accession)


def find_samples(series: list[str], library: str, species: list[str]) -> list[str]:
    def find_srr(sample: str) -> list[str]:
        print(f'Working on {sample=}')
        page = requests.get(NCBI_GEO_URL, params=dict(acc=sample))
        page.raise_for_status()

        soup = parse(page.text)

        strategy = soup.find('td', string=search_for('library strategy'))\
                       .find_parent('tr').findAll('td')[1].string
        if strategy.string.lower().strip() != library:
            return []

        organism = soup.find('td', string=search_for('organism'))\
                       .find_parent('tr').findAll('td')[1].string
        if organism.string.lower().strip() not in species:
            return []

        search = soup.find('strong', string=search_for('relations')).find_parent('tr')
        while (search := search.nextSibling) is not None:
            with contextlib.suppress(TypeError):
                if search.find('td', string=search_for('sra')):
                    sra = search.findAll('td')[1].string
                    break
        else:
            return []

        print(f'Found relevant {sra=}')

        page = requests.get(NCBI_SRA_URL, params=dict(term=sra))
        page.raise_for_status()

        return [srr.string for srr in parse(page.text).tbody.findAll('a', string=search_accession)]

    def find_samples(series: str) -> list[str]:
        print(f'Searching {series=}')
        page = requests.get(NCBI_GEO_URL, params=dict(acc=series))
        page.raise_for_status()

        soup = parse(page.text)

        samples = [
            sample.string for sample in soup.find_all(
                'a',
                href=lambda tag: tag and tag.startswith('/geo/query'),
                string=search_accession
            )
        ]
        print(f'Found {len(samples)} samples with accession numbers: {", ".join(samples)}')

        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            return list(itertools.chain.from_iterable(pool.map(find_srr, samples)))

    print(f'Searching {series=} for {library=} and {species=} samples')

    library = library.lower().strip()
    species = [ s.lower().strip() for s in species ]

    return list(itertools.chain.from_iterable(map(find_samples, series)))


def main():
    with contextlib.suppress(KeyboardInterrupt):
        import argparse

        parser = argparse.ArgumentParser(add_help=True)

        parser.add_argument('samples', nargs='+')
        parser.add_argument('-l', '--library', default='ncRNA-seq')
        parser.add_argument('-s', '--species', default=['homo_sapiens'], nargs='*', type=str)

        args = vars(parser.parse_args())
        args['species'] = [
            s.translate(str.maketrans(string.punctuation, ' ' * len(string.punctuation))).lower()
            for s in args['species']
        ]
        print(args)

        print('---BEGIN---')
        print('SEARCHING GEO FOR SERIES') 
        samples = find_samples(args.pop('samples'), args.pop('library'), args.pop('species'))

        print('DOWNLOADING APPROPRIATE SERIES')
        download_samples(samples)
        print('---END---')


if __name__ == '__main__':
    main()

