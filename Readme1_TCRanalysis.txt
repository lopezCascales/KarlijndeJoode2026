First, I tried TRUST4, a tool specifically designed to extract TCR and BCR sequences 
#from bulk RNA-seq without prior enrichment. 
# Unfortunately, TRUST4 yielded zero clonotypes across all six samples, 
# likely reflecting the limited depth of TCR-specific reads in standard bulk RNA-seq libraries, or maybe I did some mistakes.
# Then, I re-ran the TCR clonotype analysis using MiXCR (v4.7.0), 
# this time with the rna-seq preset instead of the vergani-et-al-2017-full-length preset that was used originally. 
# The original preset is designed for a specific UMI-based full-length protocol and
# is not appropriate for standard bulk RNA-seq libraries, which explains why it yielded zero clonotypes across all samples. 
# The rna-seq preset,  is specifically designed to handle fragmented, non-enriched reads, and successfully 
# recovered TCR clonotypes in all six samples (range: 10–429 per sample).


# Total clonotypes loaded: 1104
Per sample:
  # A tibble: 12 × 4
  Sample Group  Chain     N
<chr>  <fct>  <chr> <int>
  1 C01    4-1BB+ TRA      98
2 C01    4-1BB+ TRB     247
3 D01    4-1BB- TRA       5
4 D01    4-1BB- TRB      15
5 E01    4-1BB+ TRA       6
6 E01    4-1BB+ TRB       4
7 F01    4-1BB+ TRA      83
8 F01    4-1BB+ TRB     124
9 G01    4-1BB- TRA     189
10 G01    4-1BB- TRB     240
11 H01    4-1BB- TRA      35
12 H01    4-1BB- TRB      58

# The results show no significant difference in either clonotype number (Wilcoxon p=1.0) 
# or distribution between the two subsets. This supports the interpretation that the
# transcriptional differences we observe between 4-1BB+ and 4-1BB- CD4 T cells — including the IL-17 phenotype — 
# are not driven by differences in TCR clonotype composition or clonal expansion, and thus are 
# unlikely to reflect differences in antigen specificity.
