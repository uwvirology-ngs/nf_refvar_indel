process BWA_ALIGN_REF {
    tag "${meta.id}"
    label 'process_medium_java'
    container 'community.wave.seqera.io/library/bwa_samtools:eac4ad78deba8f5d'

    input:
    tuple val(meta), path(fastq)
    path ref

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"),    emit: bam
    tuple val(meta), path("*.flagstat"),                                emit: flagstat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    // Set the ref_seq variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containg a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( ref.isDirectory() ) {
        if ( ref ==~ /(.\/)?ref\/?/ ) {
            ref_seq = ''
        } else {
            ref_seq = "path=${ref}"
        }
    } else {
        ref_seq = "ref=${ref}"
    }

    """
    # index the reference genome for bwa
    bwa index ${ref}

    # generate sequence alignment map (.sam)
    bwa mem -t ${task.cpus} ${ref} ${fastq} > ${prefix}.sam

    # convert to .bam
    samtools view -bS -o ${prefix}.bam ${prefix}.sam

    # remove unmapped reads, sort and index bam files  
    samtools view -b -F 4 -@ ${task.cpus} ${prefix}.bam -o ${prefix}_mapped.bam           
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${prefix}_mapped.bam  
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam 

    # generate alignment statistics
    samtools flagstats -@ ${task.cpus} ${prefix}.sorted.bam > ${prefix}.flagstat
    """
}
