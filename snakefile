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
            --filename {output.dataset_package} --released-after 01/01/2020
        """

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
        config="config.yaml"
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
        results="results/sort_results.tsv"
    shell:
        """
        nextclade sort -m {input.minimizer} -r {output.results}  {input.data} --max-score-gap 0.5 --min-score 0.05 --min-hits 1  --all-matches
        """

rule parse_nextclade_output:
    input:
        results="results/sort_results.tsv",
        data="results/ncbi_dataset/data/genomic.fna"
    output:
        parsed=expand("results/segment_{segment}.fasta", segment=SEGMENTS)
    shell:
        """
        python scripts/parse_nextclade_sort_output.py \
            --sort-results {input.results} \
            --sequences {input.data} \
        """

# for i in {1..8}; do   
# seqkit sample -n 1000 --out-file ../results/sample_segment_${i}.fasta ../results/segment_${i}.fasta
# augur align --sequences ../results/sample_segment_${i}.fasta --output ../results/segment_${i}_aligned.fasta 
# augur tree --alignment ../results/segment_${i}_aligned.fasta --output ../results/segment_${i}.nwk 
# augur refine --tree ../results/segment_${i}.nwk --output-tree ../results/segment_${i}_refined.nwk  
# augur traits --tree ../results/segment_${i}_refined.nwk \
#             --metadata ../results/nextclade_sort.tsv \
#             --output-node-data ../results/segment${i}_traits.json \
#             --columns "inferredSegment inferredSubType ncbiSegment ncbiSubtype" \
#             --metadata-id-columns seqName 
# augur export v2 \
#             --tree ../results/segment_${i}_refined.nwk \
#             --metadata ../results/nextclade_sort.tsv \
#             --node-data ../results/segment${i}_traits.json \
#             --output ../results/segment${i}_auspice.json \
#             --metadata-id-columns "seqName"  --auspice-config auspice_config.json
# done 

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

rule all_trees:
    input:
        expand("results/segments/tree{segment}.nwk", segment=SEGMENTS)
