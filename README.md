# pcawgAdditiveVariance
>Non-coding passengers in Multiple Myeloma

This repository consist of code relevant for additive variance analysis performed on PCAWG mutations.

## Package Update 
pcawgAdditiveVariance can now fully accommodate multi-threaded execution.

All original code and documents for the pcawgAdditiveVariance suite of MATLAB scripts were developed
by [Gerstein Lab](https://github.com/gersteinlab/pcawgAdditiveVariance) and the updates
outlined in this README were incorporated by [Patrick
Blaney](https://github.com/pblaney/pcawgAdditiveVariance)

## Dependencies

Following dependencies are required to run this workflow.

* FunSeq2 (http://funseq2.gersteinlab.org/)

This is only required for the creating of `-I <funSeqOutFile>` input data of the `generateSummaryInfo.py` script, see Pre-processing step

* Python

This is required to execute the `generateSummaryInfo.py` script, specifically the library `docopt`, see Pre-processing step

* Matlab

**Required**, must have this for execution of additive variance model

* [GCTA](http://cnsgenomics.com/software/gcta/#Overview)

**Required**, the executable `gcta64` is in the `gctaFiles/` directory

## Workflow

This workflow consist of two components: Pre-processing and Post-processing

*******************

### Pre-processing step

1) Create the necessary expected directories

```
mkdir -p bedFiles
mkdir -p summaryFiles
mkdir -p SNVstats
mkdir -p results
mkdir -p keys
mkdir -p machMats
```

2) Create coding/non-coding driver gene file with the expected columns, but no column names those displayed below are for easy reading

| mutation | region | gene | cancer_abrv | cancer | chrom | start | end | strand | type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mut1 | CDS | KRAS | MM | Multiple-Myeloma | chr12 | 25205246 | 25250929 | - | snv |
| mut2 | CDS | NRAS | MM | Multiple-Myeloma | chr1 | 114704469 | 114716771 | - | snv |
| mut3 | CDS | DIS3 | MM | Multiple-Myeloma | chr13 | 72752169 | 72781900 | - | snv |


3) Create summary files for each `cohortName.null.bed` and `cohortName.null.bed` using `generateSummaryInfo.py` script. It is critical to ensure the names of samples within the two `.bed` files match, this is not checked for internally and will cause error in GCTA step

```
# Null
python generateSummaryInfo.py \
  -d codingAndNoncodingDrivers.txt \
  -I bedFiles/Multiple-Myeloma.null.bed \
  -O summaryFiles/Multiple-Myeloma.null.summary.txt

# Obs
python generateSummaryInfo.py \
  -d codingAndNoncodingDrivers.txt \
  -I bedFiles/Multiple-Myeloma.obs.bed \
  -O summaryFiles/Multiple-Myeloma.obs.summary.txt
```

*******************

### Post-processing step

1) Ensure proper input files are present

```
ls bedFiles/
cohortName.null.bed cohortName.obs.bed

ls summaryFiles/
cohortName.null.summary.txt cohortName.obs.summary.txt
```

2) Execute the pipeline
>Should be launched as a batch job using at least 6 threads 

2.1) Run get_SNVstats for null and observed separately
>Should be launched as two separate batch jobs to take advantage of parallel execution for both simultaneously

```
# Null
./get_SNVstats.sh null

# Obs
./get_SNVstats.sh obs
```

2.2) Run the remaining steps of the pipeline
>Should be launched as a batch job using at least 6 threads

```
./additive_variance.sh
```

### Results

Results text file shows calculated additive variance for each funseq threshold, which is ~0 (1e-6), along with associated p-values (0.5 indicates that no significant genetic variance was found).
Values of -1 in the results file for funseq thresholds 5 and 6 indicate that insufficient data was found at these thresholds.

### Citation
>Kumar S, Warrell J, Li S, McGillivray PD, Meyerson W, Salichos L, Harmanci A, Martinez-Fundichely A, Chan CWY, Nielsen MM, Lochovsky L, Zhang Y, Li X, Lou S, Pedersen JS, Herrmann C, Getz G, Khurana E, Gerstein MB. Passenger Mutations in More Than 2,500 Cancer Genomes: Overall Molecular Functional Impact and Consequences. Cell. 2020 Mar 5;180(5):915-927.e16. doi: 10.1016/j.cell.2020.01.032. Epub 2020 Feb 20. PMID: 32084333; PMCID: PMC7210002.
