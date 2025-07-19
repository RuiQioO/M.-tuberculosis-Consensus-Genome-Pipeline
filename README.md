# M.-tuberculosis-Consensus-Genome-Pipeline
This Nextflow pipeline maps raw paired-end FASTQ files of a single Mycobacterium tuberculosis (M. tuberculosis) sample to a reference genome, and generates a consensus genome sequence in FASTA format.

## **Features**
	- • Fully automated: from raw FASTQ to consensus FASTA
	- • Modular Nextflow DSL2 structure for clarity and scalability
	- • User can specify:
	- • Input FASTQ directory
	- • Reference genome FASTA path
	- • Output FASTA filename prefix
	- • Output directory for results
	- • Uses minimap2, samtools, and bcftools (industry-standard tools)
	- • Conda environment YAML provided for reproducibility

Data Source
This pipeline was tested with public M. tuberculosis sequencing data from NCBI SRA:
	• SRA Run Accession: SRX27529710
	• Download Example:
               fastq-dump --split-files SRX27529710

Workflow Overview
	1. Index Reference Genome
		minimap2 and samtools index the input reference FASTA.
	2. Map Reads
		Paired-end FASTQ reads are aligned to the reference using minimap2.
	3. Convert SAM to Sorted BAM
		samtools converts and sorts the alignment file.
	4. Variant Calling
		bcftools generates a VCF of detected variants.
	5. Consensus Generation
		bcftools consensus applies variants to the reference to produce a consensus genome FASTA.

Input Requirements
	• Paired-end FASTQ files in one directory
		Filenames should follow this pattern: SAMPLEID_1.fastq and SAMPLEID_2.fastq
	• Reference genome FASTA (e.g., NC_000962.3.fasta)
		Ideally, use H37Rv version 3 as specified.

Configuration
	• Set parameters in nextflow.config (example below):
		params {
   	            reads_dir  = '/path/to/fastq'
 		    ref_fasta  = '/path/to/NC_000962.3.fasta'
     		    out_prefix = 'consensus'
  	            outdir     = '/path/to/output_dir'
}

Conda Environment (Recommended)
	• To ensure full reproducibility and compatibility, use the provided Conda environment file:
		conda env create -f mtb_pipeline_env.yaml
		conda activate mtb_pipeline_env

How to Run
	1. Prepare your environment (Nextflow, minimap2, samtools, bcftools available in PATH, or use conda as above).
	2. Edit nextflow.config with your data locations.
	3. Run the pipeline:
		nextflow run main.nf
	4. Result:
		The consensus genome FASTA file will be in ${params.outdir}/consensus.fasta (or use your chosen prefix).

Output
	• consensus.fasta: The consensus genome sequence with sample-specific variants.
	• Intermediate files (.sam, .bam, .vcf.gz etc.) are managed automatically by Nextflow and stored in the work directory.

Software Requirements
	• Nextflow (version ≥ 20.10 recommended)
	• minimap2
	• samtools
	• bcftools
	• Or use the provided mtb_pipeline_env.yaml conda environment for all dependencies

Author
	• Name: Rui Qi
	• Email: rui.qi@ndm.ox.ac.uk

