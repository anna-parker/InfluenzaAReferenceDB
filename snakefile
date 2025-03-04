TAXON_ID = 11320
SEGMENTS = range(1, 9)

if os.uname().sysname == "Darwin":
    # Don't use conda-forge unzip on macOS
    # Due to https://github.com/conda-forge/unzip-feedstock/issues/16
    unzip = "/usr/bin/unzip"
else:
    unzip = "unzip"


rule fetch_ncbi_dataset_package:
    output:
        dataset_package="data/ncbi_dataset.zip",
    shell:
        """
        datasets download virus genome taxon {TAXON_ID} \
            --no-progressbar \
            --filename {output.dataset_package} \
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package=rules.fetch_ncbi_dataset_package.output.dataset_package,
    output:
        ncbi_dataset_sequences="results/ncbi_dataset/data/genomic.fna",
    params:
        unzip=unzip,
    shell:
        """
        {params.unzip} -o {input.dataset_package} -d results
        """


# Downloading the assemblies is taking too long - so I use a different approach
# rule fetch_refseq_assembly_package:
#     output:
#         dataset_package="assemblies/refseq_assembly.zip",
#     shell:
#         """
#         datasets download genome taxon {TAXON_ID} \
#             --assembly-source refseq \
#             --no-progressbar \
#             --filename {output.dataset_package} \
#         """

# rule fetch_ncbi_assembly_package:
#     output:
#         dataset_package="assemblies/genbank_assembly.zip",
#     shell:
#         """
#         datasets download genome taxon {TAXON_ID} \
#             --assembly-source genbank \
#             --filename {output.dataset_package} --released-after 01/01/2020
#         """

# rule extract_references:
#     input:
#         refseq_assembly="assemblies/refseq_assembly.zip",
#     output:
#         reference_genome="results/references.fasta",
#     shell:
#         """
#         python scripts/extract_references.py \
#             --refseq-assembly {input.refseq_assembly} \
#             --output-fasta {output.reference_genome} \
#         """

# rule extract_segments:
#     input:
#         reference_genome="results/references.fasta"
#     output:
#         fasta="results/segments/output_segment{segment}.fasta"
#     shell:
#         """
#         awk '/^>/ {{p = /_segment{wildcards.segment}$/}} p' {input.reference_genome} > {output.fasta}
#         """

rule download_references:
    input:
        config="config/config.yaml"
    output:
        reference_genome="results/references.fasta",
    shell:
        """
        python scripts/download_references.py \
            --config-file {input.config} \
            --segment-file {output.reference_genome} \
        """

rule create_minimizer:
    input:
        reference="results/references.fasta",
    output:
        minimizer="results/influenza_segments.minimizer.json"
    shell:
        """
        python scripts/create_minimizer_index.py \
            --references-fasta {input.reference} \
            --minimizer-json {output.minimizer} \
        """

rule nextclade_sort:
    input:
        minimizer="results/influenza_segments.minimizer.json",
        data="results/ncbi_dataset/data/genomic.fna"
    output:
        results="results/sort_results.tsv",
    shell:
        """
        nextclade sort -m {input.minimizer} -r {output.results}  {input.data} --max-score-gap 0.3 --min-score 0.1 --min-hits 2  --all-matches
        """

rule parse_nextclade_output:
    input:
        results="results/sort_results.tsv",
        data="results/ncbi_dataset/data/genomic.fna",
        script="scripts/parse_nextclade_sort_output.py"
    output:
        parsed=expand("results/segment_{segment}.fasta", segment=SEGMENTS),
        results="results/nextclade_sort.tsv"
    shell:
        """
        python {input.script} \
            --sort-results {input.results} \
            --sequences {input.data} \
        """

rule subsample_segments:
    input:
        fasta="results/segment_{segment}.fasta"
    output:
        subsampled_fasta="results/segments/sample_segment_{segment}.fasta"
    shell:
        """
        seqkit sample -n 5000 -2 --out-file {output.subsampled_fasta} {input.fasta}
        """

rule align_segments:
    input:
        reference_genome="results/segments/sample_segment_{segment}.fasta"
    output:
        aligned_fasta="results/segments/aligned_segment_{segment}.fasta"
    shell:
        """
        augur align --sequences {input.reference_genome} --output {output.aligned_fasta}
        """

rule segment_trees:
    input:
        reference_genome="results/segments/aligned_segment_{segment}.fasta"
    output:
        tree="results/segments/tree{segment}.nwk"
    shell:
        """
        augur tree --alignment {input.reference_genome}  \
             --output {output.tree}
        """

rule refine_trees:
    input:
        tree="results/segments/tree{segment}.nwk"
    output:
        refined_tree="results/segments/refined_tree{segment}.nwk"
    shell:
        """
        augur refine --tree {input.tree} --output-tree {output.refined_tree}
        """

rule traits:    
    input:
        tree="results/segments/refined_tree{segment}.nwk",
        metadata="results/nextclade_sort.tsv"
    output:
        traits="results/segments/traits{segment}.json"
    shell:
        """
        augur traits --tree {input.tree} --metadata {input.metadata} \
            --output-node-data {output.traits} \
            --columns "inferredSegment inferredSubType ncbiSegment ncbiSubTypeNA ncbiSubTypeHA annotation" \
            --metadata-id-columns seqName
        """

rule export:
    input:
        tree="results/segments/refined_tree{segment}.nwk",
        metadata="results/nextclade_sort.tsv",
        traits="results/segments/traits{segment}.json",
        auspice_config="config/auspice_config.json"
    output:
        auspice_json="auspice/auspice{segment}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.traits} \
            --output {output.auspice_json} \
            --metadata-id-columns "seqName"  \
            --auspice-config {input.auspice_config}
        """    

rule all_trees:
    input:
        expand("auspice/auspice{segment}.json", segment=SEGMENTS)
