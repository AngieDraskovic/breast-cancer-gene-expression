from pathlib import Path

import pandas as pd

STAR_COUNTS_DIR = Path('data/raw/star_counts')
PROCESSED_DIR = Path('data/processed')
COUNTS_MATRIX_PATH = PROCESSED_DIR / 'counts_matrix.parquet'
TPM_MATRIX_PATH = PROCESSED_DIR / 'tpm_matrix.parquet'
SAMPLE_METADATA_PATH = PROCESSED_DIR / 'sample_metadata.parquet'
EXPECTED_MATRIX_SHAPE = (60660, 1226)
TUMOR_GROUP = 'Primary Tumor'
NORMAL_GROUP = 'Solid Tissue Normal'
METASTATIC_GROUP = 'Metastatic'

def read_star_counts(file_path):
    df = pd.read_csv(
        file_path,
        sep='\t',
        comment='#'
    )

    df = df[~df['gene_id'].str.startswith('N_')]

    df = df.set_index('gene_id')
    return df[['unstranded', 'tpm_unstranded']]


def log_drop_sample(sample_id, dropped_file, kept_file, dropped_read_count, kept_read_count):
    print(f'Duplicate {sample_id}: dropping {dropped_file} ({dropped_read_count:,} reads), keeping {kept_file} ({kept_read_count:,} reads)')


def build_expression_matrices(sample_sheet):
    counts = {}
    tpm = {}
    file_names = {}
    for row in sample_sheet.itertuples():
        sample_df = read_star_counts(STAR_COUNTS_DIR / row.file_name)

        if row.sample_id in counts:
            current_read_count = counts[row.sample_id].sum()
            new_read_count = sample_df['unstranded'].sum()
            if current_read_count > new_read_count:
                log_drop_sample(row.sample_id, row.file_name, file_names[row.sample_id], new_read_count, current_read_count)
                continue
            log_drop_sample(row.sample_id, file_names[row.sample_id], row.file_name, current_read_count, new_read_count)

        counts[row.sample_id] = sample_df['unstranded']
        tpm[row.sample_id] = sample_df['tpm_unstranded']
        file_names[row.sample_id] = row.file_name

    return pd.DataFrame(counts), pd.DataFrame(tpm)

def build_sample_metadata(sample_sheet, counts_matrix):
    sample_metadata = sample_sheet.drop_duplicates(subset=['sample_id'])
    sample_metadata = sample_metadata[['sample_id', 'sample_type', 'case_id']]
    sample_metadata = sample_metadata.set_index('sample_id')
    sample_metadata = sample_metadata.loc[counts_matrix.columns]
    return sample_metadata

def sanity_checks(counts_matrix, tpm_matrix, sample_metadata):
    assert list(counts_matrix.columns) == list(tpm_matrix.columns), 'Counts and TPM matrices have different samples'
    assert list(sample_metadata.index) == list(counts_matrix.columns), 'Metadata rows are not in the same order as matrix columns'
    assert (counts_matrix >= 0).all().all(), 'Counts in count matrix are not greater or equal to 0'
    assert (tpm_matrix >= 0).all().all(), 'Values in tpm matrix are not greater or equal to 0'

    assert counts_matrix.shape == EXPECTED_MATRIX_SHAPE, 'Counts matrix does not have the expected shape'
    assert tpm_matrix.shape == EXPECTED_MATRIX_SHAPE, 'Tpm matrix does not have the expected shape'

    assert sample_metadata.groupby('sample_type').size().to_dict() == {
        NORMAL_GROUP: 113,
        TUMOR_GROUP: 1106,
        METASTATIC_GROUP: 7
    }, 'Expected group counts in sample metadata should be: Solid Tissue Normal:  113; Primary Tumor: 1106; Metastatic: 7'

    assert tpm_matrix.sum().between(999000, 1001000).all(), 'Tpm values do not sum up to ~ 1000000'

    assert sample_metadata.index.is_unique, 'Sample id should be unique in sample metadata'

    assert not tpm_matrix.isnull().any().any(), 'There should be no NaN values in tpm matrix'
    assert not counts_matrix.isnull().any().any(), 'There should be no NaN values in counts matrix'
    assert not sample_metadata.isnull().any().any(), 'There should be no NaN values in sample metadata'

    print('All sanity checks passed')

if __name__ == '__main__':
    sample_sheet = pd.read_csv('metadata/sample_sheet.tsv', sep='\t')

    if COUNTS_MATRIX_PATH.exists() and TPM_MATRIX_PATH.exists():
        counts_matrix = pd.read_parquet(COUNTS_MATRIX_PATH)   
        tpm_matrix = pd.read_parquet(TPM_MATRIX_PATH)   
    else:
        counts_matrix, tpm_matrix = build_expression_matrices(sample_sheet) 
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        counts_matrix.to_parquet(COUNTS_MATRIX_PATH)
        tpm_matrix.to_parquet(TPM_MATRIX_PATH)
        assert pd.read_parquet(COUNTS_MATRIX_PATH).equals(counts_matrix), 'Counts matrix changed when saved to Parquet'
        assert pd.read_parquet(TPM_MATRIX_PATH).equals(tpm_matrix), 'TPM matrix changed when saved to Parquet'

    sample_metadata = build_sample_metadata(sample_sheet, counts_matrix)
    sample_metadata.to_parquet(SAMPLE_METADATA_PATH)

    sanity_checks(counts_matrix, tpm_matrix, sample_metadata)