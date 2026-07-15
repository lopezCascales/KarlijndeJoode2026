C,F, G are CD137- and D, E, H CD137+
###############################################
  
  Analysis date: Wed Feb 18 14:18:39 CET 2026
Input file(s): I25-1078-D01_S83_R1_001.fastq.gz,I25-1078-D01_S83_R2_001.fastq.gz
Output file(s): D01_mixcr.alignments.vdjca
Version: 4.7.0; built=Wed Aug 07 21:19:48 CEST 2024; rev=976ba14139; lib=repseqio.v5.1
Command line arguments: align --report D01_mixcr.align.report.txt --json-report D01_mixcr.align.report.json --preset vergani-et-al-2017-full-length --save-output-file-names D01_mixcr.align.list.tsv --species hsa I25-1078-D01_S83_R1_001.fastq.gz I25-1078-D01_S83_R2_001.fastq.gz D01_mixcr.alignments.vdjca
Analysis time: 11.89s
Total sequencing reads: 1204907
Successfully aligned reads: 0 (0%)
Coverage (percent of successfully aligned):
  CDR3: 0 (NaN%)
  FR3_TO_FR4: 0 (NaN%)
  CDR2_TO_FR4: 0 (NaN%)
  FR2_TO_FR4: 0 (NaN%)
  CDR1_TO_FR4: 0 (NaN%)
  VDJRegion: 0 (NaN%)
Alignment failed: no hits (not TCR/IG?): 1201897 (99.75%)

# Alignment failed: no hits (not TCR/IG?): 1201897 (99.75%) Meaning your bulk RNA-seq reads basically have no TCR sequences that MiXCR recognizes.
Why this happens

# Standard RNA-seq libraries (polyA, rRNA-depleted, etc.) capture everything, but TCR transcripts are usually extremely rare.
# You have ~1.2 million reads per sample. For a T-cell, only a tiny fraction (~0.01–0.1%) might be actual TCR mRNA.
# That’s already borderline for MiXCR.
# The preset expects full-length V(D)J sequences, but your reads are either:
# too short to cover full VDJ
# fragmented
# not enriched for TCR (no targeted capture, no 5’RACE/UMI tagging)
# UMI stuff in the preset is probably trying to parse tags that don’t exist, so everything downstream drops to 0.


Alignment failed: absence of V hits: 187 (0.02%)
Alignment failed: absence of J hits: 2799 (0.23%)
Alignment failed: no target with both V and J alignments: 24 (0%)
Overlapped: 307950 (25.56%)
Overlapped and aligned: 0 (0%)
Overlapped and not aligned: 307950 (25.56%)
Alignment-aided overlaps, percent of overlapped and aligned: 0 (NaN%)
Realigned with forced non-floating bound: 1793914 (148.88%)
Realigned with forced non-floating right bound in left read: 30 (0%)
Realigned with forced non-floating left bound in right read: 30 (0%)
Trimming report:
  R1 reads trimmed left: 0 (0%)
  R1 reads trimmed right: 17620 (1.46%)
  Average R1 nucleotides trimmed left: 0.0
  Average R1 nucleotides trimmed right: 0.014623535260397691
  R2 reads trimmed left: 1971 (0.16%)
  R2 reads trimmed right: 3780 (0.31%)
  Average R2 nucleotides trimmed left: 0.001635810896608618
  Average R2 nucleotides trimmed right: 0.0031371715825370757
Tag parsing report:
  Execution time: 0ns
  Total reads: 1204907
  Matched reads: 1204907 (100%)
  Projection +R1 +R2: 1204907 (100%)
  For variant 0:
    For projection +R1 +R2:
      UMI:Left position: 0
      UMI:Right position: 16
      R1:Left position: 26
      R1:Right position: 51
      R2:Left position: 19
      Variants: 0
      Cost: 0
      UMI length: 16
      R1 length: 25
      R2 length: 32
======================================

library(immunarch)
Loading required package: ggplot2
Loading required package: immundata
Loading required package: dplyr

Attaching package: ‘dplyr’

The following object is masked from ‘package:Biobase’:

    combine

The following object is masked from ‘package:matrixStats’:

    count

The following objects are masked from ‘package:GenomicRanges’:

    intersect, setdiff, union

The following object is masked from ‘package:GenomeInfoDb’:

    intersect

The following objects are masked from ‘package:IRanges’:

    collapse, desc, intersect, setdiff, slice, union

The following objects are masked from ‘package:S4Vectors’:

    first, intersect, rename, setdiff, setequal, union

The following objects are masked from ‘package:BiocGenerics’:

    combine, intersect, setdiff, union

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Loading required package: duckplyr
✔ Overwriting dplyr methods with duckplyr methods.
ℹ Turn off with `duckplyr::methods_restore()`.
Loading immundata version 0.0.5

Attaching package: ‘immundata’

The following object is masked from ‘package:ggplot2’:

    annotate

Loading required package: patchwork
Hi, this is Vadim Nazarov - author of immunarch. immunarch is moving toward its
1.0 release. It will be faster, more user-friendly, and ready for publication.
Some functions will be deprecated or replaced with newer, more powerful
methods.

- Click `get_immunarch_news()` to read the latest update (what changed, what's
new, how to migrate, what's next).

- Click `list_immunarch_news()` to list all updates (latest: #1, Oct 2025).

- Migration guide: <https://immunomind.github.io/docs/tutorials/migration>

To load the package without this message:
`suppressPackageStartupMessages(library('immunarch'))`

- Vadim I. Nazarov
> data <- repLoad("D01_mixcr.clones.tsv")
# Calculate diversity metrics
repDiversity(data$clones, .method = "shannon")

== Step 1/3: loading repertoire files... ==

Processing "<initial>" ...
  -- [1/1] Parsing "D01_mixcr.clones.tsv" -- mixcr
                                                                                Rows: 0 Columns: 37                
── Column specification ────────────────────────────────────────────────────────
Delimiter: "\t"
chr (37): cloneId, readCount, readFraction, uniqueMoleculeCount, uniqueMolec...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
                                                                                  [!] Warning: zero clonotypes found, skipping

== Step 2/3: checking metadata files and merging files... ==

Processing "<initial>" ...
  -- Metadata file not found; creating a dummy metadata...

== Step 3/3: processing paired chain data... ==

Done!

Warning messages:
1: Unknown or uninitialised column: `Sample`. 
2: Unknown or uninitialised column: `Sample`. 
Error in UseMethod("select") : 
  no applicable method for 'select' applied to an object of class "NULL"
  
  
  
########################################################################################################

 java -jar mixcr.jar analyze vergani-et-al-2017-full-length   --species hsa   I25-1078-E01_S84_R1_001.fastq.gz   I25-1078-E01_S84_R2_001.fastq.gz   E01_mixcr

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report E01_mixcr.align.report.txt --json-report E01_mixcr.align.report.json --preset vergani-et-al-2017-full-length --save-output-file-names E01_mixcr.align.list.tsv --species hsa I25-1078-E01_S84_R1_001.fastq.gz I25-1078-E01_S84_R2_001.fastq.gz E01_mixcr.alignments.vdjca
The following tags and their roles will be associated with each output alignment:
  Payload tags: R1, R2
  Molecule tags: UMI(SQ)
Alignment: 0%
Alignment: 13%  ETA: 00:00:33
Alignment: 32.6%  ETA: 00:00:06
Alignment: 52.2%  ETA: 00:00:02
Alignment: 65.2%  ETA: 00:00:02
Alignment: 84.8%  ETA: 00:00:00
Alignment: 100%  ETA: 00:00:00
====================== report: align ======================
Analysis time: 12.51s
Total sequencing reads: 957401
Successfully aligned reads: 0 (0%)
Coverage (percent of successfully aligned):
  CDR3: 0 (NaN%)
  FR3_TO_FR4: 0 (NaN%)
  CDR2_TO_FR4: 0 (NaN%)
  FR2_TO_FR4: 0 (NaN%)
  CDR1_TO_FR4: 0 (NaN%)
  VDJRegion: 0 (NaN%)
Alignment failed: no hits (not TCR/IG?): 955283 (99.78%)
Alignment failed: absence of V hits: 72 (0.01%)
Alignment failed: absence of J hits: 2032 (0.21%)
Alignment failed: no target with both V and J alignments: 14 (0%)
Overlapped: 246816 (25.78%)
Overlapped and aligned: 0 (0%)
Overlapped and not aligned: 246816 (25.78%)
Alignment-aided overlaps, percent of overlapped and aligned: 0 (NaN%)
Realigned with forced non-floating bound: 1421170 (148.44%)
Realigned with forced non-floating right bound in left read: 20 (0%)
Realigned with forced non-floating left bound in right read: 20 (0%)
Trimming report:
  R1 reads trimmed left: 1 (0%)
  R1 reads trimmed right: 14050 (1.47%)
  Average R1 nucleotides trimmed left: 6.2669665061975075E-6
  Average R1 nucleotides trimmed right: 0.014676191063096863
  R2 reads trimmed left: 1636 (0.17%)
  R2 reads trimmed right: 3154 (0.33%)
  Average R2 nucleotides trimmed left: 0.0017202823059512158
  Average R2 nucleotides trimmed right: 0.00329433539342449
Tag parsing report:
  Execution time: 0ns
  Total reads: 957401
  Matched reads: 957401 (100%)
  Projection +R1 +R2: 957401 (100%)
  For variant 0:
    For projection +R1 +R2:
      UMI:Left position: 0
      UMI:Right position: 16
      R1:Left position: 26
      R1:Right position: 51
      R2:Left position: 19
      Variants: 0
      Cost: 0
      UMI length: 16
      R1 length: 25
      R2 length: 32

>>>>>>>>>>>>>>>>> mixcr refineTagsAndSort <<<<<<<<<<<<<<<<<
Running:
mixcr refineTagsAndSort --report E01_mixcr.refine.report.txt --json-report E01_mixcr.refine.report.json E01_mixcr.alignments.vdjca E01_mixcr.refined.vdjca
Correction will be applied to the following tags: UMI
Sorting will be applied to the following tags: UMI
Initialization: progress unknown
Applying correction & sorting alignments by UMI: progress unknown
Writing result: progress unknown
================ report: refineTagsAndSort ================
Analysis time: 505ms
Number of input records: 0
Number of output records: 0 (NaN%)
UMI correction report:
  UMI input groups: 0
  UMI input diversity: 0
  UMI output diversity: 0 (NaN%)
  UMI input records: 0
  UMI output records: 0 (NaN%)
  UMI mean records per tag: NaN
  UMI input core diversity: 0 (NaN%)
  UMI input core records: 0 (NaN%)
  UMI directly corrected diversity: 0 (NaN%)
  UMI directly corrected records: 0 (NaN%)
  UMI diversity filtered by tag quality: 0 (NaN%)
  UMI records filtered by tag quality: 0 (NaN%)
  UMI diversity filtered by whitelist: 0 (NaN%)
  UMI recursively corrected: 0
Filter report:
  Number of groups: 0
  Number of groups accepted: 0 (NaN%)
  Total records weight: 0
  Records weight accepted: 0 (NaN%)
  Operator #0:
    Effective threshold: -Infinity
    Nested thresholds:
      #0: 1
      #1: -Infinity

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report E01_mixcr.assemble.report.txt --json-report E01_mixcr.assemble.report.json E01_mixcr.refined.vdjca E01_mixcr.clns
No data to automatically infer minRecordsPerConsensus. Using default value: 0
Building pre-clones from tag groups: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 187ms
Final clonotype count: 0
Reads used in clonotypes, percent of total: 0 (NaN%)
Average number of reads per clonotype: NaN
Reads dropped due to the lack of a clone sequence, percent of total: 0 (NaN%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (NaN%)
Reads dropped due to low quality, percent of total: 0 (NaN%)
Reads dropped due to failed mapping, percent of total: 0 (NaN%)
Reads dropped with low quality clones, percent of total: 0 (NaN%)
Aligned reads processed: 0
Reads used in clonotypes before clustering, percent of total: 0 (NaN%)
Number of reads used as a core, percent of used: 0 (NaN%)
Mapped low quality reads, percent of used: 0 (NaN%)
Reads clustered in PCR error correction, percent of used: 0 (NaN%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (NaN%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 0
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (NaN%)
Reads dropped in post filtering: 0.0 (NaN%)
Reads filtered by tag prefix: 0 (NaN%)
Pre-clone assembler report:
  Number of input groups: 0
  Number of input groups with no assembling feature: 0
  Number of input alignments: 0
  Number of alignments with assembling feature: 0 (NaN%)
  Number of output pre-clones: 0
  Number of pre-clonotypes per group:  
  Number of assembling feature sequences in groups with zero pre-clonotypes: 0
  Number of dropped pre-clones by tag suffix conflict: 0
  Number of dropped alignments by tag suffix conflict: 0
  Number of core alignments: 0 (NaN%)
  Discarded core alignments: 0 (NaN%)
  Empirically assigned alignments: 0 (NaN%)
  Empirical assignment conflicts: 0 (NaN%)
  Tag+VJ-gene empirically assigned alignments: 0 (NaN%)
  VJ-gene empirically assigned alignments: 0 (NaN%)
  Tag empirically assigned alignments: 0 (NaN%)
  Number of ambiguous groups: 0
  Number of ambiguous tag+V/J-gene combinations: 0
  Ignored non-productive alignments: 0 (NaN%)
  Unassigned alignments: 0 (NaN%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout E01_mixcr.clns E01_mixcr.qc.txt E01_mixcr.qc.json

  Successfully aligned reads:                           0.0%   [ALERT]
  Off target (non TCR/IG) reads:                        99.78% [ALERT]
  Reads with no V or J hits:                            0.22%  [OK]
  Reads with no barcode:                                0.0%   [OK]
  Alignments that do not cover VDJRegion:               NaN    [ALERT]
  Tag groups that do not cover VDJRegion:               NaN    [ALERT]
  Barcode collisions in clonotype assembly:             NaN    [ALERT]
  Unassigned alignments in clonotype assembly:          NaN    [ALERT]
  Reads used in clonotypes:                             NaN    [ALERT]
  Alignments dropped due to low sequence quality:       NaN    [ALERT]
  Clones dropped in post-filtering:                     NaN    [ALERT]
  Alignments dropped in clones post-filtering:          NaN    [ALERT]
  Reads dropped in tags error correction and filtering: NaN    [ALERT]
  UMIs artificial diversity eliminated:                 NaN    [ALERT]
  Reads dropped in UMI error correction and whitelist:  NaN    [ALERT]
  Reads dropped in tags filtering:                      NaN    [ALERT]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones E01_mixcr.clns E01_mixcr.clones.tsv
Analysis finished successfully.


###################################################################################

 java -jar mixcr.jar analyze vergani-et-al-2017-full-length   --species hsa   I25-1078-H01_S87_R1_001.fastq.gz I25-1078-H01_S87_R2_001.fastq.gz   H01_mixcr

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report H01_mixcr.align.report.txt --json-report H01_mixcr.align.report.json --preset vergani-et-al-2017-full-length --save-output-file-names H01_mixcr.align.list.tsv --species hsa I25-1078-H01_S87_R1_001.fastq.gz I25-1078-H01_S87_R2_001.fastq.gz H01_mixcr.alignments.vdjca
The following tags and their roles will be associated with each output alignment:
  Payload tags: R1, R2
  Molecule tags: UMI(SQ)
Alignment: 0%
Alignment: 11.4%  ETA: 00:00:46
Alignment: 30.3%  ETA: 00:00:07
Alignment: 41.7%  ETA: 00:00:05
Alignment: 53%  ETA: 00:00:04
Alignment: 73.9%  ETA: 00:00:02
Alignment: 87.1%  ETA: 00:00:00
Alignment: 98.5%  ETA: 00:00:00
====================== report: align ======================
Analysis time: 15.11s
Total sequencing reads: 1658740
Successfully aligned reads: 1 (0%)
Coverage (percent of successfully aligned):
  CDR3: 1 (100%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 1653237 (99.67%)
Alignment failed: absence of V hits: 300 (0.02%)
Alignment failed: absence of J hits: 5145 (0.31%)
Alignment failed: no target with both V and J alignments: 57 (0%)
Overlapped: 450418 (27.15%)
Overlapped and aligned: 1 (0%)
Overlapped and not aligned: 450417 (27.15%)
Alignment-aided overlaps, percent of overlapped and aligned: 0 (0%)
Realigned with forced non-floating bound: 2416644 (145.69%)
Realigned with forced non-floating right bound in left read: 98 (0.01%)
Realigned with forced non-floating left bound in right read: 98 (0.01%)
TRA chains: 1 (100%)
TRA non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 0 (0%)
  R1 reads trimmed right: 24328 (1.47%)
  Average R1 nucleotides trimmed left: 0.0
  Average R1 nucleotides trimmed right: 0.014666554131449173
  R2 reads trimmed left: 3133 (0.19%)
  R2 reads trimmed right: 5150 (0.31%)
  Average R2 nucleotides trimmed left: 0.0018887830521962454
  Average R2 nucleotides trimmed right: 0.003104766268372379
Tag parsing report:
  Execution time: 0ns
  Total reads: 1658740
  Matched reads: 1658740 (100%)
  Projection +R1 +R2: 1658740 (100%)
  For variant 0:
    For projection +R1 +R2:
      UMI:Left position: 0
      UMI:Right position: 16
      R1:Left position: 26
      R1:Right position: 51
      R2:Left position: 19
      Variants: 0
      Cost: 0
      UMI length: 16
      R1 length: 25
      R2 length: 32

>>>>>>>>>>>>>>>>> mixcr refineTagsAndSort <<<<<<<<<<<<<<<<<
Running:
mixcr refineTagsAndSort --report H01_mixcr.refine.report.txt --json-report H01_mixcr.refine.report.json H01_mixcr.alignments.vdjca H01_mixcr.refined.vdjca
Correction will be applied to the following tags: UMI
Sorting will be applied to the following tags: UMI
Initialization: progress unknown
Applying correction & sorting alignments by UMI: 0%
Writing result: 0%
================ report: refineTagsAndSort ================
Analysis time: 474ms
Number of input records: 1
Number of output records: 1 (100%)
UMI correction report:
  UMI input diversity: 1
  UMI output diversity: 1 (100%)
  UMI input records: 1
  UMI output records: 1 (100%)
  UMI mean records per tag: 1
  UMI input core diversity: 1 (100%)
  UMI input core records: 1 (100%)
  UMI directly corrected diversity: 0 (0%)
  UMI directly corrected records: 0 (0%)
  UMI diversity filtered by tag quality: 0 (0%)
  UMI records filtered by tag quality: 0 (0%)
  UMI diversity filtered by whitelist: 0 (0%)
  UMI recursively corrected: 0
Filter report:
  Number of groups: 1
  Number of groups accepted: 1 (100%)
  Total records weight: 1
  Records weight accepted: 1 (100%)
  Operator #0:
    Effective threshold: 1.0
    Nested thresholds:
      #0: 1
      #1: 1

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report H01_mixcr.assemble.report.txt --json-report H01_mixcr.assemble.report.json H01_mixcr.refined.vdjca H01_mixcr.clns
No data to automatically infer minRecordsPerConsensus. Using default value: 0
Building pre-clones from tag groups: 100%
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 133ms
Final clonotype count: 0
Reads used in clonotypes, percent of total: 0 (0%)
Average number of reads per clonotype: NaN
Reads dropped due to the lack of a clone sequence, percent of total: 1 (100%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 0 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 0
Reads used in clonotypes before clustering, percent of total: 0 (0%)
Number of reads used as a core, percent of used: 0 (NaN%)
Mapped low quality reads, percent of used: 0 (NaN%)
Reads clustered in PCR error correction, percent of used: 0 (NaN%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (NaN%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 0
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (NaN%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
Pre-clone assembler report:
  Number of input groups: 1
  Number of input groups with no assembling feature: 1
  Number of input alignments: 1
  Number of alignments with assembling feature: 0 (0%)
  Number of output pre-clones: 0
  Number of pre-clonotypes per group:  
  Number of assembling feature sequences in groups with zero pre-clonotypes: 0
  Number of dropped pre-clones by tag suffix conflict: 0
  Number of dropped alignments by tag suffix conflict: 0
  Number of core alignments: 0 (0%)
  Discarded core alignments: 0 (NaN%)
  Empirically assigned alignments: 0 (0%)
  Empirical assignment conflicts: 0 (0%)
  Tag+VJ-gene empirically assigned alignments: 0 (0%)
  VJ-gene empirically assigned alignments: 0 (0%)
  Tag empirically assigned alignments: 0 (0%)
  Number of ambiguous groups: 0
  Number of ambiguous tag+V/J-gene combinations: 0
  Ignored non-productive alignments: 0 (0%)
  Unassigned alignments: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout H01_mixcr.clns H01_mixcr.qc.txt H01_mixcr.qc.json

  Successfully aligned reads:                           0.0%   [ALERT]
  Off target (non TCR/IG) reads:                        99.67% [ALERT]
  Reads with no V or J hits:                            0.33%  [OK]
  Reads with no barcode:                                0.0%   [OK]
  Alignments that do not cover VDJRegion:               100.0% [ALERT]
  Tag groups that do not cover VDJRegion:               100.0% [ALERT]
  Barcode collisions in clonotype assembly:             NaN    [ALERT]
  Unassigned alignments in clonotype assembly:          0.0%   [OK]
  Reads used in clonotypes:                             0.0%   [ALERT]
  Alignments dropped due to low sequence quality:       0.0%   [OK]
  Clones dropped in post-filtering:                     NaN    [ALERT]
  Alignments dropped in clones post-filtering:          NaN    [ALERT]
  Reads dropped in tags error correction and filtering: 0.0%   [OK]
  UMIs artificial diversity eliminated:                 0.0%   [OK]
  Reads dropped in UMI error correction and whitelist:  0.0%   [OK]
  Reads dropped in tags filtering:                      0.0%   [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones H01_mixcr.clns H01_mixcr.clones.tsv
Analysis finished successfully.

########################################################################################
head(counts_long)
# A tibble: 6 × 3
  GENEID          Sample    Counts
  <chr>           <chr>      <int>
1 ENSG00000000003 F01counts     43
2 ENSG00000000003 C01counts     85
3 ENSG00000000003 E01counts      3
4 ENSG00000000003 G01counts    396
5 ENSG00000000003 H01counts      0
6 ENSG00000000003 D01counts     85

head(counts_long)
# A tibble: 6 × 4
  GENEID          Sample    Counts Group 
  <chr>           <chr>      <int> <chr> 
1 ENSG00000000003 F01counts     43 4-1BB+
2 ENSG00000000003 C01counts     85 4-1BB+
3 ENSG00000000003 E01counts      3 4-1BB+
4 ENSG00000000003 G01counts    396 4-1BB-
5 ENSG00000000003 H01counts      0 4-1BB-
6 ENSG00000000003 D01counts     85 4-1BB-



head(counts_summary)
# A tibble: 6 × 3
  GENEID          Group  TotalCounts
  <chr>           <chr>        <int>
1 ENSG00000000003 4-1BB+         131
2 ENSG00000000003 4-1BB-         481
3 ENSG00000000005 4-1BB+           0
4 ENSG00000000005 4-1BB-           0
5 ENSG00000000419 4-1BB+         626
6 ENSG00000000419 4-1BB-         407

clone_counts <- counts_summary %>%
  group_by(Group) %>%
  summarise(NumClones = sum(TotalCounts > 0))

clone_counts

clone_counts
# A tibble: 2 × 2
  Group  NumClones
  <chr>      <int>
1 4-1BB+     20649
2 4-1BB-     21407


library(vegan)

diversity_summary <- counts_summary %>%
  group_by(Group) %>%
  summarise(Shannon = diversity(TotalCounts, index = "shannon"))

diversity_summary
Loading required package: permute
This is vegan 2.7-2
# A tibble: 2 × 2
  Group  Shannon
  <chr>    <dbl>
1 4-1BB+    7.82
2 4-1BB-    8.17


counts_summary <- counts_summary %>%
  group_by(Group) %>%
  mutate(Freq = TotalCounts / sum(TotalCounts)) %>%
  ungroup() %>%
  mutate(Category = case_when(
    Freq <= 1e-05 & Freq > 0          ~ "rare",
    Freq <= 1e-04 & Freq > 1e-05      ~ "low",
    Freq <= 0.001 & Freq > 1e-04      ~ "medium",
    Freq <= 0.01 & Freq > 0.001       ~ "high",
    Freq <= 1 & Freq > 0.01           ~ "hyperexpanded",
    TRUE                               ~ "none"
  ))

table(counts_summary$Group, counts_summary$Category)
        
          high hyperexpanded   low medium  none  rare
  4-1BB-    88             2  6039   2276 41303 13002
  4-1BB+    96             3  5750   2049 42061 12751

library(ggplot2)

ggplot(counts_summary, aes(x=Group, fill=Category)) +
  geom_bar() +
  scale_fill_brewer(palette="Set2") +
  theme_minimal() +
  labs(title="Distribución de clonotypes por grupo",
       y="Number of clones", x="Group")
    ggsave("clonotypes.svg")   
       
       
library(dplyr)
library(ggplot2)

# Seleccionar los top N clonotypes por grupo
topN <- 20

top_clones <- counts_summary %>%
  group_by(Group) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = topN) %>%
  ungroup()

# Gráfico de barras
ggplot(top_clones, aes(x=reorder(GENEID, Freq), y=Freq, fill=Group)) +
  geom_bar(stat="identity") +
  facet_wrap(~Group, scales="free_x") +
  coord_flip() +  # barras horizontales
  theme_minimal() +
  labs(title = paste("Top", topN, "clonotypes per group"),
       x = "Clonotype (GENEID)", y = "Relative Frequency") +
  scale_fill_manual(values=c("4-1BB+"="#1f77b4", "4-1BB-"="#ff7f0e"))       
  
  ggsave("top_clones.svg")
  
  ggplot(top_clones, aes(x=reorder(GENEID, Freq), y=Freq, fill=Group)) +
  geom_bar(stat="identity") +
  facet_wrap(~Group, scales="fixed") +  # mismo eje X para ambos
  coord_flip() +                        # barras horizontales
  theme_minimal() +
  labs(title = paste("Top", topN, "clonotypes por grupo"),
       x = "Clonotype (GENEID)", y = "Frecuencia relativa") +
  scale_fill_manual(values=c("4-1BB+"="#1f77b4", "4-1BB-"="#ff7f0e"))
  
    ggsave("top_clones_sameXscale.svg")
    
    
    library(biomaRt)
library(dplyr)

# Conectar al Ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Obtener mapeo ENSG → gene symbol
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = top_clones$GENEID,
  mart = ensembl
)

# Unir con tus datos
top_clones <- top_clones %>%
  left_join(gene_map, by = c("GENEID" = "ensembl_gene_id")) %>%
  mutate(GeneLabel = ifelse(hgnc_symbol == "", GENEID, hgnc_symbol))

# Luego usar GeneLabel en el plot
ggplot(top_clones, aes(x=reorder(GeneLabel, Freq), y=Freq, fill=Group)) +
  geom_bar(stat="identity") +
  facet_wrap(~Group, scales="fixed") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste("Top", topN, "clonotypes por grupo"),
       x = "Clonotype", y = "Relative Frequency") +
  scale_fill_manual(values=c("4-1BB+"="#1f77b4", "4-1BB-"="#ff7f0e"))
  
     ggsave("top_clones_sameXscale_symbol.svg")
     
     topN <- 50

top_clones <- counts_summary %>%
  group_by(Group) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = topN) %>%
  ungroup()
  
      library(biomaRt)
library(dplyr)

# Conectar al Ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Obtener mapeo ENSG → gene symbol
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = top_clones$GENEID,
  mart = ensembl
)

# Unir con tus datos
top_clones <- top_clones %>%
  left_join(gene_map, by = c("GENEID" = "ensembl_gene_id")) %>%
  mutate(GeneLabel = ifelse(hgnc_symbol == "", GENEID, hgnc_symbol))

# Luego usar GeneLabel en el plot
ggplot(top_clones, aes(x=reorder(GeneLabel, Freq), y=Freq, fill=Group)) +
  geom_bar(stat="identity") +
  facet_wrap(~Group, scales="fixed") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste("Top", topN, "clonotypes por grupo"),
       x = "Clonotype", y = "Relative Frequency") +
  scale_fill_manual(values=c("4-1BB+"="#1f77b4", "4-1BB-"="#ff7f0e"))
  
     ggsave("top_clones_sameXscale_symbol_top50.svg")
     
     
  
