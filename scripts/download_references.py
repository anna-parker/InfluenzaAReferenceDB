from Bio import Entrez, SeqIO
from dataclasses import dataclass
import click
import yaml

# Set your email address (required by NCBI)
Entrez.email = "your_email@example.com"

@dataclass
class SubTypes:
    references: dict[str, str]



def download_references(subtypes, segment_file):

    # Fetch the sequence in FASTA format
    with open(segment_file, "w", encoding="utf-8") as output_file:
        for name, accession in subtypes.references.items():
            with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                record = SeqIO.read(handle, "fasta")
                output_file.write(f">{name}\n{record.seq}\n")


@click.command()
@click.option("--config-file", required=True, type=click.Path(exists=True))
@click.option("--segment-file", required=False, type=click.Path())
def main(config_file: str, segment_file) -> None:

    with open(config_file, encoding="utf-8") as file:
        full_config = yaml.safe_load(file)
        relevant_config = {key: full_config[key] for key in SubTypes.__annotations__}
        subtypes = SubTypes(**relevant_config)
    download_references(subtypes, segment_file)


if __name__ == "__main__":
    main()
    