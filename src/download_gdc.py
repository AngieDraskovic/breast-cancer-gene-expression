"""Download TCGA-BRCA STAR-Counts files and their sample sheet from the GDC API.

Usage (from the project root):
    python src/download_gdc.py --metadata-only   # sample sheet only
    python src/download_gdc.py                   # sample sheet + all count files
"""

import argparse
import hashlib
import time
from pathlib import Path

import pandas as pd
import requests

FILES_ENDPOINT = 'https://api.gdc.cancer.gov/files'
DATA_ENDPOINT = 'https://api.gdc.cancer.gov/data'

RAW_DIR = Path('data/raw')
COUNTS_DIR = RAW_DIR / 'star_counts'
SAMPLE_SHEET_PATH = Path('metadata/sample_sheet.tsv')
MAX_ATTEMPTS = 5

FILTERS = {
    'op': 'and',
    'content': [
        {'op': 'in', 'content': {'field': 'cases.project.project_id', 'value': ['TCGA-BRCA']}},
        {'op': 'in', 'content': {'field': 'data_type', 'value': ['Gene Expression Quantification']}},
        {'op': 'in', 'content': {'field': 'analysis.workflow_type', 'value': ['STAR - Counts']}},
        {'op': 'in', 'content': {'field': 'access', 'value': ['open']}},
    ],
}

FIELDS = [
    'file_id',
    'file_name',
    'file_size',
    'md5sum',
    'cases.submitter_id',
    'cases.samples.submitter_id',
    'cases.samples.sample_type',
]


def fetch_sample_sheet():
    params = {
        'filters': FILTERS,
        'fields': ','.join(FIELDS),
        'format': 'JSON',
        'size': 5000,
    }
    response = requests.post(FILES_ENDPOINT, json=params, timeout=60)
    response.raise_for_status()
    hits = response.json()['data']['hits']

    rows = []
    for hit in hits:
        # Each STAR-Counts file comes from exactly one case and one sample.
        assert len(hit['cases']) == 1 and len(hit['cases'][0]['samples']) == 1, hit['file_id']
        case = hit['cases'][0]
        sample = case['samples'][0]
        rows.append({
            'file_id': hit['file_id'],
            'file_name': hit['file_name'],
            'case_id': case['submitter_id'],
            'sample_id': sample['submitter_id'],
            'sample_type': sample['sample_type'],
            'file_size': hit['file_size'],
            'md5sum': hit['md5sum'],
        })

    return pd.DataFrame(rows).sort_values('sample_id').reset_index(drop=True)


def md5_of(path):
    digest = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(row):
    """Download one file; skip it if it's already there and intact."""
    path = COUNTS_DIR / row.file_name
    if path.exists() and md5_of(path) == row.md5sum:
        return 'skipped'

    # The GDC server occasionally drops connections, so retry with a growing pause.
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = requests.get(f'{DATA_ENDPOINT}/{row.file_id}', timeout=120)
            response.raise_for_status()
            break
        except requests.exceptions.RequestException:
            if attempt == MAX_ATTEMPTS:
                raise
            time.sleep(10 * attempt)
    path.write_bytes(response.content)

    if md5_of(path) != row.md5sum:
        path.unlink()
        raise ValueError(f'Checksum mismatch for {row.file_name}')
    return 'downloaded'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata-only', action='store_true')
    args = parser.parse_args()

    SAMPLE_SHEET_PATH.parent.mkdir(parents=True, exist_ok=True)
    sample_sheet = fetch_sample_sheet()
    sample_sheet.to_csv(SAMPLE_SHEET_PATH, sep='\t', index=False)

    print(f'Saved {SAMPLE_SHEET_PATH} ({len(sample_sheet)} files)')
    print(sample_sheet['sample_type'].value_counts().to_string())
    print(f'Unique patients: {sample_sheet['case_id'].nunique()}')
    print(f'Total size: {sample_sheet['file_size'].sum() / 1e9:.1f} GB')

    if args.metadata_only:
        return

    COUNTS_DIR.mkdir(parents=True, exist_ok=True)
    for i, row in enumerate(sample_sheet.itertuples(), start=1):
        status = download_file(row)
        print(f'[{i}/{len(sample_sheet)}] {status}: {row.file_name}')


if __name__ == '__main__':
    main()
