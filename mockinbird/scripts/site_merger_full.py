import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('all_sites')
    parser.add_argument('--factor_table')
    parser.add_argument('mock_table')
    parser.add_argument('out_table')
    args = parser.parse_args()

    df_ref = pd.read_table(args.all_sites)

    if args.factor_table:
        df_factor = pd.read_table(args.factor_table)

        combined = df_ref.merge(df_factor, on=['chrom', 'pos', 'strand'], how='left')
        del df_factor
        df_mock = pd.read_table(args.mock_table)
        combined = combined.merge(df_mock, on=['chrom', 'pos', 'strand'], how='left',
                                  suffixes=('_factor', '_mock'))
        del df_mock
        col_names = [
            'chrom', 'pos',
            'k_factor', 'n_factor',
            'k_mock', 'n_mock',
            'strand',
        ]
        int_cols = [
            'k_factor', 'n_factor',
            'k_mock', 'n_mock',
        ]
    else:
        df_mock = pd.read_table(args.mock_table)
        df_mock = df_mock.rename(columns={'k': 'k_mock', 'n': 'n_mock'})
        combined = df_ref.merge(df_mock, on=['chrom', 'pos', 'strand'], how='left')
        col_names = [
            'chrom', 'pos',
            'k_mock', 'n_mock',
            'strand',
        ]
        int_cols = [
            'k_mock', 'n_mock',
        ]

    final = combined[col_names]

    final.fillna(0, inplace=True)
    final[int_cols] = final[int_cols].astype(int)

    final.to_csv(args.out_table, index=False, sep='\t')


if __name__ == '__main__':
    main()
