## Intro ##

This is a pipeline for a full variant calling work-flow starting with fastq files to be downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1172971<br>

These fastq files come from raw sequencing reads from a breast cancer cell line called JIMT1. The purpose of this pipeline is to identify variants between a reference and JIMT1, and find genes that are the most heavily mutated. I make use of a gene annotation file (gtf) called `gencode.gtf.gz` which can be found in this folder. 

## To Run ##

IMPORTANT: `pipeline.sh` will probably take around 24-48 hours to run, and is not necessary to see the results of the analysis. It is more meant to document how I was able to get the output.vcf.gz and all the QC and alignment steps up to it.<br>

In order to run the `variants.py` file, which contains the post-variant calling analysis, do the following:

- Install dependencies from requirements.txt with `pip install -r requirements.txt`
- Make sure you have python 3.12 or greater (older version might work too, but this is best for reproducibility)
- Run `python variants.py`

This will run all of the analysis that is seen in Report.pdf
