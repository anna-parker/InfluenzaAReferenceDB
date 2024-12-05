import logging
import os
import re

import click
from Bio import SeqIO
import zipfile


logger = logging.getLogger(__name__)
logging.basicConfig(
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)8s (%(filename)20s:%(lineno)4d) - %(message)s ",
    datefmt="%H:%M:%S",
)

def extract_subtype(seq_name):
    pattern = r'\bH([1-9]|1[0-8])N([1-9]|1[0-2])\b'  # Match H1-18 and N1-12
    match = re.search(pattern, seq_name)
    return match.group(0) if match else ""

def extract_segment(seq_name):
    pattern = r'segment ([1-8])'
    match = re.search(pattern, seq_name)
    return match.group(0) if match else ""


def search_assembly_dir(assembly_dir, subtype_reference_map, subtypes):
    for dirpath, dirnames, filenames in os.walk(assembly_dir):
        for filename in filenames:
            if not filename.endswith('.fna'):
                continue
            file_path = os.path.join(dirpath, filename)
            print(f'File: {file_path}')
            annotated_subtypes = set()
            annotated_segments = {}
            with open(file_path, encoding="utf-8") as f_in:
                records = SeqIO.parse(f_in, "fasta")
                for record in records:
                    annotated_subtypes.add(extract_subtype(record.description))
                    annotated_segments[extract_segment(record.description)] = record.seq
            if len(annotated_subtypes) == 1 and annotated_subtypes <= subtypes and len(annotated_segments) == 8:
                subtype = annotated_subtypes.pop()
                for segment in annotated_segments:
                    subtype_reference_map[(subtype, segment)] = annotated_segments[segment]
                subtypes.remove(subtype)
                print("Found a reference for subtype: ", subtype)
    return subtype_reference_map, subtypes


@click.command(help="Parse fasta header, only keep if fits regex filter_fasta_headers")
@click.option("--refseq-assembly-path", required=True, type=click.Path(exists=True))
@click.option("--genbank-assembly-path", required=False, type=click.Path(exists=True))
@click.option("--output-fasta", required=True, type=str)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
)
def main(
    refseq_assembly_path: str,
    genbank_assembly_path: str,
    output_fasta: str,
    log_level: str,
) -> None:
    logger.setLevel(log_level)

    subtypes = set(["H"+ str(i) + "N" + str(j) for i in range(1, 19) for j in range(1, 12)])
    subtype_reference_map = {}
    print(f"Searching for references for {len(subtypes)} subtypes")

    with zipfile.ZipFile(refseq_assembly_path, 'r') as zip_ref:
        zip_ref.extractall(refseq_assembly_path.split('.')[0])
    refseq_assembly_dir = os.path.join(os.curdir, refseq_assembly_path.split('.')[0])
    subtype_reference_map, subtypes = search_assembly_dir(refseq_assembly_dir, subtype_reference_map, subtypes)
    print(subtype_reference_map)

    print(f"Still searching for references for {len(subtypes)} subtypes")

    with zipfile.ZipFile(genbank_assembly_path, 'r') as zip_ref:
        zip_ref.extractall(refseq_assembly_path.split('.')[0])
    refseq_assembly_dir = os.path.join(os.curdir, refseq_assembly_path.split('.')[0])
    subtype_reference_map, subtypes = search_assembly_dir(refseq_assembly_dir, subtype_reference_map, subtypes)
    
    print(f"Still searching for references for {len(subtypes)} subtypes")

    with open(output_fasta, "w", encoding="utf-8") as output_file:
        for subtype, segment in subtype_reference_map.items():
            name = f"{subtype[0]}_{subtype[1]}".replace(" ", "")
            output_file.write(f">{name}\n{segment}\n")

    


if __name__ == "__main__":
    main()