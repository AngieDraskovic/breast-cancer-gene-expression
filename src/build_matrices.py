import pandas as pd


STAR_COUNTS_DIR = 'data/raw/star_counts/'


def read_star_counts(file_path):
    df = pd.read_csv(
        file_path,
        sep='\t',
        comment='#'
    )

    df = df[~df['gene_id'].str.startswith('N_')]

    df = df.set_index('gene_id')
    return df['unstranded']


def log_drop_sample(sample_id, dropped_file, kept_file, dropped_read_count, kept_read_count):
    print(f'Duplicate {sample_id}: dropping {dropped_file} ({dropped_read_count:,} reads), keeping {kept_file} ({kept_read_count:,} reads)')


def build_counts_matrix(sample_sheet):
    counts = {}
    file_names = {}
    for row in sample_sheet.itertuples():
        new_series = read_star_counts(STAR_COUNTS_DIR + row.file_name)
        if row.sample_id in counts:
            current_read_sum = sum(counts[row.sample_id])
            new_read_sum = sum(new_series)
            if current_read_sum > new_read_sum:
                dropped_file = row.file_name
                dropped_read_count = new_read_sum
                kept_read_count = current_read_sum
                log_drop_sample(row.sample_id, dropped_file, file_names[row.sample_id], dropped_read_count, kept_read_count)
                continue
            else:
                dropped_file = file_names[row.sample_id]
                dropped_read_count = current_read_sum
                kept_read_count = new_read_sum
                log_drop_sample(row.sample_id, dropped_file, row.file_name, dropped_read_count, kept_read_count)

        counts[row.sample_id] = new_series
        file_names[row.sample_id] = row.file_name

    matrix = pd.DataFrame.from_dict(counts)

    return matrix


if __name__ == '__main__':
    sample_sheet = pd.read_csv(
        'metadata/sample_sheet.tsv',
        sep='\t',
    )
    print(build_counts_matrix(sample_sheet))