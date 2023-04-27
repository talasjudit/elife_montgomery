Calling SNPs:
Input - FASTA of reference genome with only autosomes, indexed with bowtie2 and samtools, .bam or .fastq.gz files of whole-genome sequencing of the alternate genome.
1. dna.snps.sh
2. snpsplit.prep.sh
Output - FASTA of reference genome with only autosomes with all SNP positions as Ns, .txt file of SNPs for use with SNPSplit.

Mapping RNA-seq data
Input - .bam files of each RNA-seq replicate, tab-delimited .txt file of input bam files and sample names, indexed genome as done when calling SNPs.
1. mapping.sh
Output - Mapped .bam, expression counts per gene .genes.results.

Chromatin profiling
Input - Indexed genome, indexed spike-in genome, merged .bam files of each CUT&RUN replicate, tab-delimited .txt file of input bam files and sample names, .bed file of Protein Coding Gene coordinates.
1. cutrun.sh
2. bedtools_map.sh
Output - Mapped .bam, normalized .bw coverage, chromatin enrichment per gene.

Assignment of reads to parental genome of origin
Input - Mapped .bam files
1. snpsplit.run.sh
Output - .bam files of reads sorted based on parental genome of origin

FISH-IF analyses
Input - .tif files of separate channels of IF and FISH and DAPI
1. sexchrFIS.ijm
2. H3K27me3FISH.ijm OR H3K9me1FISH.ijm
3. FISH.R
Output - .csv tables of overlaps between foci, plots of percentage overlap between immunostain foci and FISH signal foci.

Nuclei IF analyses
Input - 
1. deconvolution using Huygens and settings in .hgsd file
2. Milos
3. NucleiIF.R
Output - Deconvolved images, .csv table of measurements of signal and mask overlaps, plots of these measurements.

Measurement of wild type and mutant development
Input - .csv and .xls .txt tables of measurements of WT and mutant plants
1. gemmae_area.ijm
2. embryo_measurements.R
Output - Plots of comparisons between WT and mutant plants.

Generate most plots in R
Input - Output tables from above steps, TPM values from .genes.results files as single table.
1. OmniPrep.v3.R
2. OmniPlotv3.R
Output - All the plots you'll ever need ;)