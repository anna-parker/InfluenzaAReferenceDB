### Annotating Influenza A

For Genspectrum we would like to be able to annotate all influenza A sequences with their segment and their subtype.

This is not trivial as NCBI no longer assigns new influenza sequences to a subtype taxon and NCBI annotations in file names and qualifiers are not standardized and sometimes have errors. 

Here I take the approach suggested by @cornelius-roemer and use [nextclade sort](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/reference.html#nextclade-sort) to perform fast local alignment of sequence k-mers to find the highest scoring sequence matches. To use nextclade sort I need to create a minimizer-index, this must contain references for all the sequences that I would like to annotate. In this case this means I need at least a reference for each segment, but also for all the subtypes. In influenza A subtypes are determined by the HA and NA segments, there are a total of 18 HA subtypes and 11 NA subtypes. To annotate the HA and NA subtypes I use reference sequences used previously in the literature. 

For HA I relied solely on the excellently annotated:
__Abdulrahman DA, Meng X, Veit M. S-Acylation of Proteins of Coronavirus and Influenza Virus: Conservation of Acylation Sites in Animal Viruses and DHHC Acyltransferases in Their Animal Reservoirs. Pathogens. 2021 May 29;10(6):669.__ [link](https://pmc.ncbi.nlm.nih.gov/articles/PMC8227752/#app1-pathogens-10-00669)

For NA I used a combination of
__Wohlbold TJ, Krammer F. In the shadow of hemagglutinin: a growing interest in influenza viral neuraminidase and its role as a vaccine antigen. Viruses. 2014 Jun 23;6(6):2465-94__ [link](https://pubmed.ncbi.nlm.nih.gov/24960271/) and __Jang YH, Seong BL. The Quest for a Truly Universal Influenza Vaccine. Front Cell Infect Microbiol. 2019 Oct 10;9:344.__ [link](https://pubmed.ncbi.nlm.nih.gov/31649895/). Wohlbold et al. supply a long list of NA reference segments and I use the visualization aid in Jang et al.'s paper to condense the set. Specifically 
I use only two of the listed N2 references, choosing two N2 references from each subcluster found by Jang et al., and I chose a H5N1-N1 reference that is the center of the N1 cluster. 

I then create auspice trees using a subsample of the annotated segments for a visual verification. I additionally plot the NCBI assigned reference if it was contained in the sequence description in a parsable format. From visual inspection my annotation appears to be significantly better than the available NCBI annotations for the NA segment and comparable for the HA segment.