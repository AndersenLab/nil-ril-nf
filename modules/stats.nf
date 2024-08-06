
/*
    fq bam stats
*/

process fq_bam_stats {

    tag { ID }

    label "trim"
    label "xs"

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")

    output:
        file 'bam_stat'

    """
        cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_bam_stats {

    container null
    executor "local"

    publishDir params.out + "/fq", mode: 'copy'

    input:
        val stat_files

    output:
        file("fq_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
        cat ${stat_files.join(" ")} >> fq_bam_stats.tsv
    """
}

/*
    SM bam stats
*/

process SM_bam_stats {

    label "trim"

    tag { SM }

    input:
        tuple file(SM), file("${SM}.bam"), file("${SM}.bam.bai")

    output:
        file 'bam_stat' 

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    container null
    executor "local"

    publishDir params.out + "/SM", mode: 'copy'

    input:
        val stat_files

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}

/*
 SM idx stats
*/

process idx_stats_SM {
    
    label "trim"

    tag { SM }

    input:
        tuple file(SM), file("${SM}.bam"), file("${SM}.bam.bai")
    output:
        file 'SM_bam_idxstats'

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > SM_bam_idxstats
    """
}

process combine_idx_stats {

    container null
    executor "local"

    publishDir params.out +"/SM", mode: 'copy'

    input:
        val bam_idxstats

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
    """

}
