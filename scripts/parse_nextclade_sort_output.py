import re
import pandas as pd
import numpy as np
from Bio import SeqIO
import click


def extract_subtype(seq_name):
    pattern = r"\(H([1-9]|1[0-8])N([1-9]|1[0-2])\)"  # Match H1-18 and N1-12
    match = re.search(pattern, seq_name)
    return match.group(0) if match else None

def extract_subtypeNA(seq_name):
    if extract_subtype(seq_name):
        pattern = r"N([1-9]|1[0-2])\)"
        match = re.search(pattern, extract_subtype(seq_name))
        return match.group(0)[:-1] if match else np.nan
    return np.nan

def extract_subtypeHA(seq_name):
    if extract_subtype(seq_name):
        pattern = r"H([1-9]|1[0-8])N"
        match = re.search(pattern, extract_subtype(seq_name))
        return match.group(0)[:-1] if match else np.nan
    return np.nan


def extract_segment(seq_name):
    pattern = r"segment ([1-8])"
    match = re.search(pattern, seq_name)
    return match.group(0).replace("ment ", "") if match else np.nan


def extract_accession(seq_name):
    return seq_name.split(" ")[0]


def extract_computed_subtype(dataset):
    if extract_computed_segment(dataset) in ["seg4", "seg6"]:
        return dataset.split("_")[1] if len(dataset.split("_")) > 1 else np.nan
    return np.nan

def extract_computed_annotation(dataset):
    if extract_computed_segment(dataset) not in ["seg4", "seg6"]:
        return dataset.split("_")[2] if len(dataset.split("_")) > 2 else np.nan
    return np.nan


def extract_computed_segment(dataset):
    return dataset.split("_")[0]


def parse_file(sort_results, sequences):
    df = pd.read_csv(sort_results, sep="\t", dtype={"index": "Int64"})

    no_rows = df.shape[0]

    # Ensure 'score' is treated as a numeric column (replace 'score' with the actual column name if different)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    # Drop rows where 'score' is NaN (optional, if there are invalid scores)
    df = df.dropna(subset=["score"])

    print(f"Dropped {no_rows - df.shape[0]} rows with NaN score")

    # Group by 'index', then sort within each group by 'score' and keep the highest score
    df_sorted = df.sort_values(["index", "score"], ascending=[True, False])
    df_highest_per_group = df_sorted.drop_duplicates(subset="index", keep="first")

    df_highest_per_group.loc[:, "ncbiSubTypeNA"] = df_highest_per_group["seqName"].apply(
        extract_subtypeNA
    )
    df_highest_per_group.loc[:, "ncbiSubTypeHA"] = df_highest_per_group["seqName"].apply(
        extract_subtypeHA
    )
    df_highest_per_group.loc[:, "ncbiSegment"] = df_highest_per_group["seqName"].apply(
        extract_segment
    )
    df_highest_per_group.loc[:, "seqName"] = df_highest_per_group["seqName"].apply(
        extract_accession
    )
    df_highest_per_group.loc[:, "inferredSubType"] = df_highest_per_group["dataset"].apply(
        extract_computed_subtype
    )
    df_highest_per_group.loc[:, "inferredSegment"] = df_highest_per_group["dataset"].apply(
        extract_computed_segment
    )
    df_highest_per_group.loc[:, "annotation"] = df_highest_per_group["dataset"].apply(
        extract_computed_annotation
    )
    df_highest_per_group.drop(columns="dataset", inplace=True)

    # Reset index and save the result (optional)
    df_highest_per_group.reset_index(drop=True, inplace=True)

    # Optionally save to a new file
    output_path = "results/nextclade_sort.tsv"
    df_highest_per_group.to_csv(output_path, sep="\t", index=False)
    non_null_count = df_highest_per_group["ncbiSegment"].notna().sum()
    print(f"Number of non-null rows in column 'ncbiSegment': {non_null_count}")

    filtered_df = df_highest_per_group[
        (df_highest_per_group["inferredSegment"] != df_highest_per_group["ncbiSegment"])
        & (df_highest_per_group["ncbiSegment"].notna())
    ]

    # Count the number of such rows
    count_mismatch = len(filtered_df)

    # Display the count
    print(
        f"Number of rows where 'inferredSegment' is not equal to 'ncbiSegment' and 'ncbiSegment' is not null: {count_mismatch}"
    )
    print(filtered_df)

    segment_lookup = {}
    for _, row in df_highest_per_group.iterrows():
        segment_lookup[row["seqName"]] = row["inferredSegment"]

    # Create file handles for each segment
    output_files = {
        segment: open(f"results/segment_{segment}.fasta", "w", encoding="utf-8")
        for segment in range(1, 9)
    }

    try:
        # Parse the genomic file once
        with open(sequences, encoding="utf-8") as f_in:
            for record in SeqIO.parse(f_in, "fasta"):
                inferred_segment = segment_lookup.get(record.id)
                if inferred_segment:
                    segment_number = int(inferred_segment[-1])
                    if 1 <= segment_number <= 8:
                        output_files[segment_number].write(
                            f">{record.id}\n{record.seq}\n"
                        )
    finally:
        # Ensure all files are closed properly
        for file in output_files.values():
            file.close()


@click.command()
@click.option("--sort-results", required=True, type=click.Path(exists=True))
@click.option("--sequences", required=False, type=click.Path(exists=True))
def main(sort_results: str, sequences) -> None:
    parse_file(sort_results, sequences)


if __name__ == "__main__":
    main()
