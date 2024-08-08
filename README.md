# pcawgAdditiveVariance
>Extension of PCAWG additive variance analysis of non-coding passengers to Multiple Myeloma genomes

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

* Python (Successful with v3.9.13)

This is required to execute the `generateSummaryInfo.py` script, specifically the library `docopt`, see Pre-processing step

* MATLAB (Successful with R2023)

**Required**, must have this for execution of additive variance model

* [GCTA](http://cnsgenomics.com/software/gcta/#Overview)

**Required**, the executable `gcta64` is in the `gctaFiles/` directory

## Workflow

This workflow consist of two components: Pre-processing and Post-processing

*******************

### Pre-processing step

1) Create the necessary expected directories

```
mkdir -p bedFiles summaryFiles SNVstats results keys machMats
```

2) Create coding/non-coding driver gene file, the `generateSummaryInfo.py` script expects the format displayed below, however final driver gene file should not have a header

| mutation | region | gene | cancer_abrv | cancer | chrom | start | end | strand | type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mut1 | CDS | KRAS | MM | Multiple-Myeloma | chr12 | 25205246 | 25250929 | - | snv |
| mut2 | CDS | NRAS | MM | Multiple-Myeloma | chr1 | 114704469 | 114716771 | - | snv |
| mut3 | CDS | DIS3 | MM | Multiple-Myeloma | chr13 | 72752169 | 72781900 | - | snv |


3) Create summary files for each `cohortName.null.bed` and `cohortName.null.bed` using `generateSummaryInfo.py` script. It is critical to ensure the names of samples within the two `.bed` files match, this is not checked for internally and will cause error in GCTA step

For the **null** dataset summary file
```
python generateSummaryInfo.py \
  -d codingAndNoncodingDrivers.txt \
  -I bedFiles/Multiple-Myeloma.null.bed \
  -O summaryFiles/Multiple-Myeloma.null.summary.txt
```

For the **obs** dataset summary file
```
python generateSummaryInfo.py \
  -d codingAndNoncodingDrivers.txt \
  -I bedFiles/Multiple-Myeloma.obs.bed \
  -O summaryFiles/Multiple-Myeloma.obs.summary.txt
```

*******************

### Post-processing step

1) Ensure proper input files are present

The FunSeq2 BED files
```
ls bedFiles/
```
>Example: Multiple-Myeloma.null.bed Multiple-Myeloma.obs.bed

The null and obs summary files
```
ls summaryFiles/
```
>Example: Multiple-Myeloma.null.summary.txt Multiple-Myeloma.obs.summary.txt

2) Execute the pipeline as set of multi-threaded batch jobs

Run get_SNVstats for null and observed separately, roughly ~2hr per set
For **null**
```
./get_SNVstats.sh null
```

<details>
<summary>SLURM Example</summary>
<br>

```
sbatch \
  --job-name=getSNVstats \
  --time=4:00:00 \
  --mem=24G \
  --cpus-per-task=6 \
  --wrap="./get_SNVstats.sh null"
```

</details>

For **obs**
```
./get_SNVstats.sh obs
```

<details>
<summary>SLURM Example</summary>
<br>

```
sbatch \
  --job-name=getSNVstats \
  --time=4:00:00 \
  --mem=24G \
  --cpus-per-task=6 \
  --wrap="./get_SNVstats.sh obs"
```

</details>

Run the remaining steps of the pipeline, may take up to 24hrs
```
./additive_variance.sh
```

<details>
<summary>SLURM Example</summary>
<br>

```
sbatch \
  --job-name=addVar \
  --time=36:00:00 \
  --mem=36G \
  --cpus-per-task=6 \
  --wrap="./additive_variance.sh"
```

</details>

### Results

Results text file shows calculated additive variance for each funseq threshold, which is ~0 (1e-6), along with associated p-values (0.5 indicates that no significant genetic variance was found).
Values of -1 in the results file for funseq thresholds 5 and 6 indicate that insufficient data was found at these thresholds.

### Citation
>Kumar S, Warrell J, Li S, McGillivray PD, Meyerson W, Salichos L, Harmanci A, Martinez-Fundichely A, Chan CWY, Nielsen MM, Lochovsky L, Zhang Y, Li X, Lou S, Pedersen JS, Herrmann C, Getz G, Khurana E, Gerstein MB. Passenger Mutations in More Than 2,500 Cancer Genomes: Overall Molecular Functional Impact and Consequences. Cell. 2020 Mar 5;180(5):915-927.e16. doi: 10.1016/j.cell.2020.01.032. Epub 2020 Feb 20. PMID: 32084333; PMCID: PMC7210002.
