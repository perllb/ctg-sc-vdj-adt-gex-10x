# ctg-sc-vdj-adt-gex-10x
Nextflow pipeine for 10x VDJ+ADT+GEX libraries - based on cellranger

- Supports mm10 and hg38 references, but can also be run with custom reference genome and annotation (must be added via nextflow.config). See custom genome below.

## USAGE

1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-rna-10x/tree/master/container/sc-rna-10x.v6
2. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
3. Edit the nextflow.config file to fit your project and system. 
4. Run pipeline 
```
nohup nextflow run pipe-sc-vdj-adt-gex-10x.nf > log.pipe-sc-vdj-adt-gex-10x.txt &
```

## Input Files

The following files must be in the runfolder to start pipeline successfully.

1. Samplesheet (`CTG_SampleSheet.sc-vdj-adt-gex-10x.csv`)

### Samplesheet requirements:

1. Above row with [Data], declare metaid (name/id of the current analysis), email (for customer delivery) and featureref (currently only `totalseqc` is supported).
2. The "gex" samples defines the sample ids. other libraries (tcr, bcr and fb) must have the same sample_ID, but added "_TCR", "_BCR" or "_FB" suffix. 
3. Sample_Pair groups the "samples" that belong together across GEX, TCR and BCR.
4. Leave Lane column empty if not applied, otherwise set to 1,2 etc.
 
Note: One samplesheet pr project!
Note: Must be in comma-separated values format (.csv)

| , | , | , | , | , | , | , |
| --- | --- | --- | --- | --- | --- | --- |
| metaid | 2021_999 | , | , | , | , | , |
| email | customer@email.com | , | , | , | , | , |
| featureref | totalseqc | , | , | , | , | , |
| [Data] | , | , | , | , | , | , |
| **Lane** | **Sample_ID** | **index** | **Sample_Project** | **Sample_Species** | **Sample_Lib** | **Sample_Pair** |
|  | Si1 | SI-GA-D9 | 2021_012 | human | gex | 1 |
|  | Si2 | SI-GA-H9 | 2021_012 | human | gex | 2 |
|  | Si1_TCR | SI-GA-G9 | 2021_012 | human | tcr | 1 |
|  | Si2_TCR | SI-GA-I9 | 2021_012 | human | tcr | 2 |
|  | Si1_BCR | SI-GA-J9 | 2021_012 | human | bcr | 1 |
|  | Si2_BCR | SI-GA-K9 | 2021_012 | human | bcr | 2 |
|  | Si1_FB | SI-GA-L9 | 2021_012 | human | fb | 1 |
|  | Si2_FB | SI-GA-M9 | 2021_012 | human | fb | 2 |


The nf-pipeline takes the following Columns from samplesheet to use in channels:

- `Sample_ID` : ID of sample. Sample_ID can only contain a-z, A-Z and "_".  E.g space and hyphen ("-") are not allowed! If 'Sample_Name' is present, it will be ignored. The "gex" samples defines the sample ids. other libraries (tcr, bcr and fb) must have the same sample_ID, but added "_TCR", "_BCR" or "_FB" suffix. 
- `index` : Must use index ID (10x ID) if dual index. For single index, the index sequence works too.
- `Sample_Project` : Project ID. E.g. 2021_033, 2021_192_userA.
- `Sample_Species` : Only 'human'/'mouse'/'custom' are accepted. This will set reference genome and annotation gtf, as well as vdj reference. If species is not human or mouse, set 'custom'. This custom reference genome has to be specified in the nextflow config file. See below how to edit the config file.
- `Lane` : Lane of flowcell for this sample. Leave Lane column empty if not applied, otherwise set to 1,2 etc. 

The metadata above [Data] is used by driver.
- `metaid` is optional
- `featureref`: only totalseqc is supported.

### Samplesheet template (.csv)

#### Name : `CTG_SampleSheet.sc-vdj-adt-gex-10x.csv`
```
metaid,2021_081
email,anna.hagstrom@med.lu.se
featureref,totalseqc
[Data]
Lane,Sample_ID,Sample_Name,index,Sample_Project,Sample_Species,Sample_Lib,Sample_Pair
,45D_viable_cells,21_147,SI-TT-E8,2021_081,human,gex,1
,113D_Leukemia_cells,21_148,SI-TT-F8,2021_081,human,gex,2
,113D_CD45pos_CD19neg,21_149,SI-TT-G8,2021_081,human,gex,3
,45D_viable_cells_BCR,21_147_BCR,SI-TT-B8,2021_081,human,bcr,1
,113D_Leukemia_cells_BCR,21_148_BCR,SI-TT-C8,2021_081,human,bcr,2
,113D_CD45pos_CD19neg_BCR,21_149_BCR,SI-TT-D8,2021_081,human,bcr,3
,45D_viable_cells_FB,21_147_FB,SI-TN-A2,2021_081,human,fb,1
,113D_Leukemia_cells_FB,21_148_FB,SI-TN-B2,2021_081,human,fb,2
,113D_CD45pos_CD19neg_FB,21_149_FB,SI-TN-C2,2021_081,human,fb,3
``` 


## Pipeline steps:

Cellranger version: cellranger v6.0 

* `Demultiplexing` (cellranger mkfastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/mkfastq).
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `Align` + `Counts` (cellranger multi): Aligns fastq files to reference genome, counts genes for each cell/barcode, performs VDJ and Feature Barcode analysis combined with GEX, secondary analysis such as clustering and generates the cloupe files (https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi).
* `Aggregation` (cellranger aggr): Automatically creates the input csv pointing to molecule_info.h5 files for each sample to be aggregated and executes aggregation (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate). This is only run if there is more than one sample pr project. **NOT YET SUPPORTED**
* `Cellranger count metrics` (bin/ctg-sc-vdj-adt-gex-count-metrics-concat.py): Collects main count metrics (#cells and #reads/cell etc.) from each sample and collect in table
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report (https://multiqc.info/)
 * `md5sum`: md5sum of all generated files


## Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * cellranger metrics: Main metrics summarising the count / cell output 
        * fastqc output 
        * multiqc output: Summarizing Cellranger, FastQC and demultiplexing 
    * `fastq`: Contains raw fastq files from cellranger mkfastq.
    * `count`: Cellranger multi output. Here you find gene/cell count matrices for GEX, VDJ and Feature Barcoding, secondary analysis output, and more. See (https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi) for more information on the output files.
    * `summaries`: 
        * web-summary files which provide an overview of essential metrics from the 10x run. 
        * cloupe files which can be used to explore the data interactively in the Loupe browser (https://support.10xgenomics.com/single-cell-vdj/software/visualization/latest/what-is-loupe-vdj-browser)  
    * `aggregate`:
        * Output from cellranger aggregation. This is only run if there is more than one sample pr project. **NOT YET SUPPORTED**
    * `ctg-md5.PROJ_ID.txt`: text file with md5sum recursively from output dir root    


## Container
https://github.com/perllb/ctg-containers/tree/main/sc-rna-10x/sc-rna-10x.v6

## Custom genome 

If custom genome (not hg38 or mm10) is used

1. Set "Sample_Species" column to 'custom' in samplesheet for samples of custom references :

Example:
| , | , | , | , | , | , | , |
| --- | --- | --- | --- | --- | --- | --- |
| ProjectID | 2021_999 | , | , | , | , | , |
| email | customer@email.com | , | , | , | , | , |
| featureref | totalseqc | , | , | , | , | , |
| [Data] | , | , | , | , | , | , |
| **Lane** | **Sample_ID** | **index** | **Sample_Project** | **Sample_Species** | **Sample_Lib** | **Sample_Pair** |
|  | Si1 | SI-GA-D9 | 2021_012 | **custom** | gex | 1 |
|  | Si2 | SI-GA-H9 | 2021_012 | human | gex | 2 |
|  | Si1_TCR | SI-GA-G9 | 2021_012 | **custom** | tcr | 1 |
|  | Si2_TCR | SI-GA-I9 | 2021_012 | human | tcr | 2 |
|  | Si1_BCR | SI-GA-J9 | 2021_012 | **custom** | bcr | 1 |
|  | Si2_BCR | SI-GA-K9 | 2021_012 | human | bcr | 2 |
|  | Si1_FB | SI-GA-L9 | 2021_012 | **custom** | fb | 1 |
|  | Si2_FB | SI-GA-M9 | 2021_012 | human | fb | 2 |
 
 2. In nextflow.config, set 
 `custom_genome=/PATH/TO/CUSTOMGENOME`
 
## Add custom genes (e.g. reporters) to cellranger annotation

You can use this script to add custom genes to the cellranger ref
https://github.com/perllb/ctg-cellranger-add2ref
