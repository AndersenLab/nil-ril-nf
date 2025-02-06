
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
    label "sm"

    publishDir params.out + "/bam", mode: 'copy'

    tag { SM }

    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
        tuple val(SM), path(bams), path(indices)

    output:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), emit: merged_SM
        path "${SM}.duplicates.txt", emit: duplicates_file

    """
    count=`echo ${bams.join(" ")} | tr ' ' '\\n' | wc -l`
    if [ "\${count}" -eq "1" ]; then
        mv ${bams[0]} ${SM}.merged.bam
        mv ${indices[0]} ${SM}.merged.bam.bai
    else
        samtools merge -@ ${task.cpus} ${SM}.merged.bam ${bams.join(" ")}
        samtools index -@ ${task.cpus} ${SM}.merged.bam
    fi
    picard MarkDuplicates I=${SM}.merged.bam \\
                            O=${SM}.bam \\
                            M=${SM}.duplicates.txt \\
                            VALIDATION_STRINGENCY=LENIENT \\
                            REMOVE_DUPLICATES=true
    samtools index -@ ${task.cpus} ${SM}.bam
    """
}

