
/*
    Fastq coverage
*/
process fq_coverage {

    tag { ID }

    label "bam"
    label "xs"

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")
    output:
        file("${ID}.coverage.tsv")


    """
        bam coverage ${ID}.bam > ${ID}.coverage.tsv
    """
}

process fq_coverage_merge {

    container null
    executor "local"

    publishDir params.out + "/fq", mode: 'copy'

    input:
        path fq_set

    output:
        file("fq_coverage.full.tsv")
        path "fq_coverage.tsv", emit: fq_coverage_plot

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> fq_coverage.full.tsv
        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}

/*
    Coverage Bam
*/
process SM_coverage {

    label "bam"
    label "xs"

    tag { SM }

    input:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai") 

    output:
        file("${SM}.coverage.tsv")

    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process SM_coverage_merge {

    container null
    executor "local"

    publishDir params.out + "/SM", mode: 'copy'

    input:
        val sm_set 

    output:
        file("SM_coverage.full.tsv")
        path "SM_coverage.tsv", emit: SM_coverage_plot

    """
        echo -e 'SM\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv
        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
    """

}
