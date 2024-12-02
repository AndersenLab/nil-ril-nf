
process generate_issue_plots {

    label "R"
    label "xs"

    publishDir params.out + "/plots", mode: 'copy'

    errorStrategy 'ignore'

    input:
        tuple file("fq_coverage.tsv"), file("SM_bam_idxstats.tsv"), file("SM_coverage.tsv"), file("bam_duplicates.tsv")

    output:
        file("coverage_comparison.png")
//        file("coverage_comparison.svg")
        file("unmapped_reads.png")
//        file("unmapped_reads.svg")
        file("duplicates.png")
        file("duplicates.pdf")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/qc_plots.R
    """
}

process output_tsv {

    label "annotation"
    label "sm"

    publishDir params.out + "/hmm", mode: 'copy'

    input:
        tuple file("NIL.hmm.vcf.gz"), file("NIL.hmm.vcf.gz.csi")

    output:
        file("gt_hmm_genotypes.tsv")

    """
        cat <(echo -e "CHROM\tPOS\tSAMPLE\tGT\tGT_ORIG\tAD") <(bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\t%GT_ORIG\t%AD\n]' NIL.hmm.vcf.gz | sed 's/0\\/0/0/g' | sed 's/1\\/1/1/g') > gt_hmm_genotypes.tsv
    """

}

// generate unique breakpoint sites for RIL cross object
process generate_cross_object {

    label "vcfkit"
    label "sm"

    publishDir params.out + "/cross_object", mode: 'copy'

    input:
        tuple file("NIL.hmm.vcf.gz"), file("NIL.hmm.vcf.gz.csi"), file("SM_coverage.tsv"), file("gt_hmm_fill.tsv")

    output:
        //file("cross_obj.Rdata")
        //file("comparegeno.png")
        //file("comparegeno.Rda")
        //file("rug.png")
        //file("estrf.Rda")
        //file("CM_map.Rda")
        //file("identicals_list.Rda")
        file("cross_obj_geno.tsv")
        //file("cross_obj_pheno.tsv")
        file("cross_obj_strains.tsv")
        file("breakpoint_sites.tsv.gz")

    """
    # Generate breakpoint sites
    cat <(cut -f 1,2 gt_hmm_fill.tsv) <(cut -f 1,3 gt_hmm_fill.tsv) |\
    grep -v 'chrom' |\
    sort -k1,1 -k2,2n |\
    uniq > breakpoint_sites.tsv
    bgzip breakpoint_sites.tsv -c > breakpoint_sites.tsv.gz && tabix -s1 -b2 -e2 breakpoint_sites.tsv.gz
    # Generate output strains list
    awk  '\$2 != "coverage" { print }'  SM_coverage.tsv  |\
    cut -f 1 |\
    grep -v '${params.A}' |\
    grep -v '${params.B}' |\
    sort -V > cross_obj_strains.tsv
    paste <(echo -e "marker\\tchrom\\tpos") <(cat cross_obj_strains.tsv| tr '\n' '\t' | sed 's/\t\$//g') > cross_obj_geno.tsv
    bcftools view -T breakpoint_sites.tsv.gz -m 2 -M 2 NIL.hmm.vcf.gz |\
    vk filter MISSING --max=0.00001 - |\
    bcftools query --samples-file cross_obj_strains.tsv -f '%CHROM\\_%POS\\t%CHROM\\t%POS[\\t%GT]\n' |\
    awk  -v OFS='\t' '''
            {   
                gsub("0/0", "0", \$0);
                gsub("1/1", "1", \$0);
                gsub("./.","", \$0);
                \$3;
                print
            }
        ''' - >> cross_obj_geno.tsv
    
    # Rscript --vanilla `which generate_cross_object.R`
    """
}
