nextflow.enable.dsl=2

params.reads_dir  = params.reads_dir ?: "./fastq"
params.ref_fasta  = params.ref_fasta ?: "./ref/NC_000962.3.fasta"
params.out_prefix = params.out_prefix ?: "consensus"
params.outdir     = params.outdir ?: "./result"

// 1. 自动配对fastq，支持 _R1/_R2 命名
Channel
    .fromFilePairs("${params.reads_dir}/*_{1,2}.fastq", flat: true)
    .set { reads_ch }
reads_ch.view { "READS DEBUG: ${it}, type: ${it.getClass()}" }
// 2. ref fasta channel
Channel
    .fromPath(params.ref_fasta)
    .ifEmpty { error "Reference fasta not found: ${params.ref_fasta}" }
    .set { ref_fasta_ch }

// 3. 对参考基因组建索引
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

// 4. 比对
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

// 5. sam 转 sorted bam
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

// 6. 变异检测
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

// 7. consensus 生成
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

// ----------- workflow 主体部分 -----------
workflow {

    // 参考基因组索引
    index_ref_out = index_ref(ref_fasta_ch)
    mmi_ch = index_ref_out.mmi

    // 生成 mapping 输入：[sample_id, reads], ref_fasta, mmi
    mapping_input_ch = reads_ch
        .combine(ref_fasta_ch)
        .combine(mmi_ch)
        .map { it.flatten() } // [sample_id, reads, ref_fasta, mmi]

    mapped_ch = mapping(mapping_input_ch)

    // sam2sortedbam 输入： [sample_id, sam]
    bam_ch = sam2sortedbam(mapped_ch)

    // variant_calling 输入：[sample_id, sorted_bam], ref_fasta
    variant_input_ch = bam_ch
        .combine(ref_fasta_ch)
        .map { it.flatten() } // [sample_id, sorted_bam, ref_fasta]

    vcf_ch = variant_calling(variant_input_ch)

    // consensus 输入：[sample_id, vcf_gz], ref_fasta
    consensus_input_ch = vcf_ch
        .combine(ref_fasta_ch)
        .map { it.flatten() }

    consensus(consensus_input_ch)
 }

