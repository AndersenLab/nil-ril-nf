
process format_duplicates {

    container null
    executor "local"

    publishDir params.out + "/duplicates", mode: 'copy'

    input:
        val duplicates_set 

    output:
        file("bam_duplicates.tsv")

    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}

process merge_bam {

    label "trim"
    label "custom"

    tag { SM }

    errorStrategy 'retry'
    time { 2.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    publishDir params.out + "/bam", mode: 'copy'

    input:
        tuple val(SM), path(bams), path(indices)

    output:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), emit: merged

    """
    samtools merge -@ ${task.cpus} ${SM}.bam ${bams.join(" ")}
    samtools index -@ ${task.cpus} ${SM}.bam
    """
}

process mark_duplicates {

    label "alignment"
    label "custom"

    errorStrategy 'retry'
    time { 2.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }
    //publishDir params.out + "/bam", mode: 'copy'

    tag { SM }

    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
        tuple val(SM), path("duplicate_${SM}.bam"), path(index)

    output:
        tuple val(SM), file("${SM}.bam"), emit: deduped
        path "${SM}.duplicates.txt", emit: duplicates_file

    """
    picard MarkDuplicates I=duplicate_${SM}.bam \\
                            O=${SM}.bam \\
                            M=${SM}.duplicates.txt \\
                            VALIDATION_STRINGENCY=SILENT \\
                            REMOVE_DUPLICATES=true \\
                            TAGGING_POLICY=All \\
                            REMOVE_SEQUENCING_DUPLICATES=true
    """
}


process index_bam {

    label "trim"
    label "xs"

    //publishDir params.out + "/bam", mode: 'copy'

    tag { SM }

    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
        tuple val(SM), path("${SM}.bam")

    output:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), emit: bam

    """
    samtools index -@ ${task.cpus} ${SM}.bam
    """
}
