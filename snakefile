TAXON_ID = 11320
SEGMENTS = range(1, 9)


rule fetch_refseq_assembly_package:
    output:
        dataset_package="assemblies/refseq_assembly.zip",
    shell:
        """
        datasets download genome taxon {TAXON_ID} \
            --assembly-source refseq \
            --no-progressbar \
            --filename {output.dataset_package} \
        """

rule fetch_ncbi_assembly_package:
    output:
        dataset_package="assemblies/genbank_assembly.zip",
    shell:
        """
        datasets download genome taxon {TAXON_ID} \
            --assembly-source genbank \
            --no-progressbar \
            --filename {output.dataset_package} \
        """

rule extract_references:
    input:
        refseq_assembly="assemblies/refseq_assembly.zip",
        genbank_assembly="assemblies/genbank_assembly.zip",
    output:
        reference_genome="results/references.fasta",
    shell:
        """
        python scripts/extract_references.py \
            --refseq-assembly {input.refseq_assembly} \
            --genbank-assembly {input.genbank_assembly} \
            --output-fasta {output.reference_genome} \
        """

# As annotations are sometimes wrong, visualize segments by aligning and building trees
rule extract_segments:
    input:
        reference_genome="results/references.fasta"
    output:
        fasta="results/segments/output_segment{segment}.fasta"
    shell:
        """
        awk '/^>/ {{p = /_segment{wildcards.segment}$/}} p' {input.reference_genome} > {output.fasta}
        """

rule align_segments:
    input:
        reference_genome="results/segments/output_segment{segment}.fasta"
    output:
        aligned_fasta="results/segments/aligned_segment{segment}.fasta"
    shell:
        """
        augur align --sequences {input.reference_genome} --output {output.aligned_fasta}
        """

rule segment_trees:
    input:
        reference_genome="results/segments/aligned_segment{segment}.fasta"
    output:
        tree="results/segments/tree{segment}.nwk"
    shell:
        """
        augur tree --alignment {input.reference_genome}  \
             --output {output.tree}
        """


# nextclade sort -m influenza_segments.minimizer.json -r out.tsv  sample_influenza.fasta --max-score-gap 0.5 --min-score 0.05 --min-hits 1  --all-matches