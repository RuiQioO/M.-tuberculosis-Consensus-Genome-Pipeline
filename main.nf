// Enable Nextflow DSL2
nextflow.enable.dsl=2

// ====== PARAMETERS ======
params.reads_dir  = params.reads_dir ?: "./fastq"
params.ref_fasta  = params.ref_fasta ?: "./ref/NC_000962.3.fasta"
params.out_prefix = params.out_prefix ?: "consensus"
params.outdir     = params.outdir ?: "./result"

// ====== 1. CREATE READS CHANNEL ======
// Pair R1/R2 fastq files for a single sample (naming: *_1.fastq, *_2.fastq)
Channel
    .fromFilePairs("${params.reads_dir}/*_{1,2}.fastq", flat: true)
    .set { reads_ch }
reads_ch.view { "READS DEBUG: ${it}, type: ${it.getClass()}" }

// ====== 2. CREATE REFERENCE FASTA CHANNEL ======
Channel
    .fromPath(params.ref_fasta)
    .ifEmpty { error "Reference fasta not found: ${params.ref_fasta}" }
    .set { ref_fasta_ch }

// ====== 3. INDEX REFERENCE GENOME ======
// This process generates minimap2 index (.mmi) and samtools fai index (.fai) for the reference
process index_ref {
    tag "indexing"
    input:
        path ref_fasta
    output:
        path("${ref_fasta}.mmi"), emit: mmi
        path("${ref_fasta}.fai"), emit: fai
    script:
    """
    minimap2 -d ${ref_fasta}.mmi ${ref_fasta}
    samtools faidx ${ref_fasta}
    """
}

// ====== 4. MAP READS TO REFERENCE ======
// Aligns paired-end reads to reference using minimap2
process mapping {
    tag "${sample_id}"
    input:
        tuple val(sample_id), path(reads1),path(reads2), path(ref_fasta), path(mmi)
    output:
        tuple val(sample_id), path("${sample_id}.sam")
    script:
    """
    mkdir -p ${params.outdir}
    minimap2 -ax sr ${mmi} ${reads1} ${reads2} > ${sample_id}.sam
    """
}

// ====== 5. CONVERT SAM TO SORTED BAM ======
// Converts SAM to sorted BAM, then indexes BAM file
process sam2sortedbam {
    tag "${sample_id}"
    input:
        tuple val(sample_id), path(sam)
    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam")
    script:
    """
    mkdir -p ${params.outdir}
    samtools view -Sb ${sam} | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

// ====== 6. VARIANT CALLING (BCFTOOLS) ======
// Uses bcftools to call variants (VCF + index), ready for consensus step
process variant_calling {
    tag "${sample_id}"
    input:
        tuple val(sample_id), path(sorted_bam), path(ref_fasta)
    output:
        tuple val(sample_id), path("${sample_id}.vcf.gz"),path("${sample_id}.vcf.gz.tbi")
    script:
    """
    mkdir -p ${params.outdir}
    bcftools mpileup -Ou -f ${ref_fasta} ${sorted_bam} | bcftools call -c -Oz -o ${sample_id}.vcf.gz
    bcftools index --tbi ${sample_id}.vcf.gz
    """
}

// ====== 7. CONSENSUS FASTA GENERATION ======
// Applies variants to reference to create consensus genome FASTA
process consensus {
    tag "${sample_id}"
    input:
        tuple val(sample_id), path(vcf_gz), path(tbi), path(ref_fasta)
    output:
        path("${params.out_prefix}.fasta")
    script:
    """
    mkdir -p ${params.outdir}
    cat ${ref_fasta} | bcftools consensus ${vcf_gz} > ${params.outdir}/${params.out_prefix}.fasta
    """
}

// ====== WORKFLOW DEFINITION ======
// Connects all channels and processes in order
workflow {

    // Index the reference genome (minimap2 + samtools)
    index_ref_out = index_ref(ref_fasta_ch)
    mmi_ch = index_ref_out.mmi

    // Combine fastq, reference fasta, and index for mapping
    mapping_input_ch = reads_ch
        .combine(ref_fasta_ch)
        .combine(mmi_ch)
        .map { it.flatten() } // [sample_id, reads, ref_fasta, mmi]

    mapped_ch = mapping(mapping_input_ch)

    // Convert SAM to sorted BAM
    bam_ch = sam2sortedbam(mapped_ch)

    // Variant calling: input BAM and reference
    variant_input_ch = bam_ch
        .combine(ref_fasta_ch)
        .map { it.flatten() } // [sample_id, sorted_bam, ref_fasta]

    vcf_ch = variant_calling(variant_input_ch)

    // Consensus: input VCF, its index, and reference fasta
    consensus_input_ch = vcf_ch
        .combine(ref_fasta_ch)
        .map { it.flatten() }

    // Generate consensus FASTA in the output directory
    consensus(consensus_input_ch)
 }

