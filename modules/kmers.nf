
/*
    Kmer counting
*/
process kmer_counting {

    label 'alignment'
    label "sm"

    tag { ID }

    input:
        tuple val(SM), val(ID), val(LB), file(fq1), file(fq2), path(genome_path), val(genome_basename)
    output:
        file("${ID}.kmer.tsv")

    """
    # fqs will have same number of lines
    export OFS="\t"
    fq_wc="`zcat ${fq1} | awk 'NR % 4 == 0' | wc -l`"
    zcat ${fq1} ${fq2} | fastq-kmers -k 6 | awk -v OFS="\t" -v ID=${ID} -v SM=${SM} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${ID}.kmer.tsv
    """
}

process merge_kmer {
    
    container null
    executor "local"

    publishDir params.out + "/phenotype", mode: 'copy'

    input:
        file("kmer*.tsv")

    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tSM\tID\twc") *.tsv > kmers.tsv
    """

}
