
process output_hmm {

    label "vcfkit"

    publishDir params.out + "/hmm", mode: 'copy'

    input:
        tuple file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi")

    output:
        file("gt_hmm.tsv")

    """
        vk hmm --transition=${params.transition} --A=${params.A} --B=${params.B} NIL.filter.vcf.gz > gt_hmm.tsv
    """

}

process output_hmm_fill {

    label "vcfkit"

    publishDir params.out + "/hmm", mode: 'copy'

    input:
        tuple file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi") 

    output:
        file("gt_hmm_fill.tsv") 

    """
        vk hmm --transition=${params.transition} --infill --endfill --A=${params.A} --B=${params.B} NIL.filter.vcf.gz > gt_hmm_fill.tsv
    """

}


process output_hmm_vcf {

    label "vcfkit"

    publishDir params.out + "/hmm", mode: 'copy'

    input:
        tuple file("NIL.vcf.gz"), file("NIL.vcf.gz.csi") 

    output:
        tuple file("NIL.hmm.vcf.gz"), file("NIL.hmm.vcf.gz.csi") 

    """
        vk hmm --transition=${params.transition} --vcf-out --A=${params.A} --B=${params.B} NIL.vcf.gz | bcftools view -O z > NIL.hmm.vcf.gz
        bcftools index NIL.hmm.vcf.gz
    """

}


process plot_hmm {

    label "R"

    publishDir params.out + "/hmm", mode: 'copy'

    input:
        file("gt_hmm_fill.tsv")

    output:
        file("gt_hmm.png")
        file("gt_hmm.pdf")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/plot_hmm.R "${params.cA}" "${params.cB}"
    """

}
