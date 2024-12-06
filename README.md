# Annotating Influenza A

For Genspectrum we would like to be able to annotate all influenza A sequences with their segment and their subtype.

This is not trivial as NCBI no longer assigns new influenza sequences to a subtype taxon and NCBI annotations in file names and qualifiers are not standardized and sometimes have errors. 

Here I take the approach suggested by @cornelius-roemer and use [nextclade sort](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/reference.html#nextclade-sort) to perform fast local alignment of sequence k-mers to find the highest scoring sequence matches. To use nextclade sort I need to create a [minimizer-index](https://github.com/nextstrain/nextclade_data/blob/master/scripts/lib/minimizer.py), this must contain references for all the sequences that I would like to annotate. In this case this means I need at least a reference for each segment, but also for all the subtypes. In influenza A subtypes are determined by the HA and NA segments, there are a total of 18 HA subtypes and 11 NA subtypes. To annotate the HA and NA subtypes I use reference sequences used previously in the literature. 

For HA I relied solely on the excellently annotated:
__Abdulrahman DA, Meng X, Veit M. S-Acylation of Proteins of Coronavirus and Influenza Virus: Conservation of Acylation Sites in Animal Viruses and DHHC Acyltransferases in Their Animal Reservoirs. Pathogens. 2021 May 29;10(6):669.__ [link](https://pmc.ncbi.nlm.nih.gov/articles/PMC8227752/#app1-pathogens-10-00669)

For NA I used a combination of
__Wohlbold TJ, Krammer F. In the shadow of hemagglutinin: a growing interest in influenza viral neuraminidase and its role as a vaccine antigen. Viruses. 2014 Jun 23;6(6):2465-94__ [link](https://pubmed.ncbi.nlm.nih.gov/24960271/) and __Jang YH, Seong BL. The Quest for a Truly Universal Influenza Vaccine. Front Cell Infect Microbiol. 2019 Oct 10;9:344.__ [link](https://pubmed.ncbi.nlm.nih.gov/31649895/). Wohlbold et al. supply a long list of NA reference segments and I use the visualization aid in Jang et al.'s paper to condense the set. Specifically 
I use only two of the listed N2 references, choosing two N2 references from each subcluster found by Jang et al., and I chose the H5N1-N1 reference that was the center of the N1 cluster for the N1 reference. 

I was initially uncertain which references to use for the other segments, so I first used the full reference assembly for the subtypes H5N1, H1N1, H2N2, H3N2, H7N9 and H9N2 to see the diversity of the other segments. Here I also see some visual clades, normally there are about three larger clades corresponding to the reference sequences from H5N1, H1N1 and H3N2 - which I chose to use in my minimizer index to improve hits (these trees can be viewed in the pre-work folder).

I then created auspice trees using a subsample of the annotated segments for a visual verification. I additionally annotated trees with the NCBI assigned subtype (`ncbiSubTypeHA` and `ncbiSubTypeNA`) if it was contained in the sequence description in a parsable format. From visual inspection my annotation appears to be significantly better than the available NCBI annotations for the NA segment and HA segment.

### Running the analysis

You can rerun the analysis using:
```
micromamba create -f environment.yml
micromamba activate influenza-db
snakemake all_trees
```

Alternatively you can run `auspice view` in the auspice directory or drop the individual auspice.json files into https://auspice.us/ to visualize the results. 