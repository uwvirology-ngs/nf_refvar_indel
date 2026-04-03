process GATK_INDELREALIGNER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'
    // gatk=3.8, samtools=1.23.1

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta5), path(known_vcf)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known = known_vcf ? "-known ${known_vcf}" : ""

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK IndelRealigner] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    samtools faidx "${fasta}"
    samtools dict "${fasta}" > "${fasta.simpleName}.dict"

    gatk3 \\
        -Xmx${avail_mem}M \\
        -T IndelRealigner \\
        -R ${fasta} \\
        -I ${bam} \\
        --targetIntervals ${intervals} \\
        ${known} \\
        -o ${prefix}.bam \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK IndelRealigner] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """
}
