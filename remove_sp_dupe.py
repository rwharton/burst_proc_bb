import pandas as pd
import click

@click.command()
@click.option("--csv_file", type=str, 
              help="Full path of csv file", required=True)
def remove_dupes(csv_file):
    df = pd.read_csv(csv_file, delim_whitespace=True) 

    rows = []
    sample_column = 3
    snr_column = 1
    downfact_column = 4
    df = df.sort_values(df.columns[sample_column])
    i = 0
    while i < len(df) - 1:
        sample_current = df.iloc[i][sample_column]
        snr_current = df.iloc[i][snr_column]
        downfact_current = df.iloc[i][downfact_column]
        max_snr_row = df.iloc[i]

        # next sample
        j = i+1

        while j < len(df) and abs(sample_current - df.iloc[j][sample_column]) <= downfact_current:
            if df.iloc[j][snr_column] > snr_current:
                max_snr_row = df.iloc[j]
                snr_current = df.iloc[j][snr_column]

            j += 1

        rows.append(max_snr_row)
        i = j

    # Check last sample
    if df.iloc[-1]['Sample'] not in pd.concat(rows)['Sample'].values:
        rows.append(df.iloc[-1])
    
    df_new = pd.concat(rows, axis=1).transpose()
    df_new = df_new.sort_values(df_new.columns[snr_column], ascending=False)
    filename_parts = csv_file.split('.')
    new_filename = f"{filename_parts[0]}_nodupes.{filename_parts[1]}"
    df_new.to_csv(new_filename, sep=' ', index=False)

debug = 0
if __name__ == "__main__":
    if not debug: 
        remove_dupes()
