
/* 
    Call variants at the individual level for concordance
*/

process split_fq {

    label "vcfkit"
    label "sm"

    tag { SM }

    input:
        tuple val(SM), val(ID), val(LIB), file("input.bam"), file("input.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi'), path(genome_path), val(genome_basename)

    output:
        tuple file('rg_gt.tsv'), val(SM)

    """
        INDEX=`find -L ${genome_path} -name "${genome_basename}.amb" | sed 's/\\.amb\$//'`
        # Split bam file into individual read groups; Ignore MtDNA
        contigs="`samtools view -H input.bam | awk '\$0 ~ "^@SQ" { gsub("SN:", "", \$2); print \$2 }' | grep -v 'MtDNA'`"
        rg_list=("`samtools view -H input.bam | awk '\$0 ~ "^@RG"' | grep -o -w -e "ID:[A-Z|a-z|0-9|_|-]*" | sed 's/ID://g'`")
        samtools split -@ ${task.cpus} -M 100 -f '%!.%.' input.bam
        # DO NOT INDEX ORIGINAL BAM; ELIMINATES CACHE!
        bam_list="`ls -1 *.bam | grep -v 'input.bam'`"
        ls -1 *.bam | grep -v 'input.bam' | xargs -I {} -P ${task.cpus} sh -c "samtools index {}"
        # Call a union set of variants
        for rg in \$rg_list; do
            echo \${contigs} | tr ' ' '\\n' | \
                xargs -I {} -P ${task.cpus} sh -c "bcftools mpileup --redo-BAQ -r {} -O b --annotate DP,AD,ADF,ADR,SP --fasta-ref \${INDEX} \${rg}.bam | \
                bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z > {}.\${rg}.vcf.gz"
            order=`echo \${contigs} | tr ' ' '\\n' | awk -v rg=\${rg} '{ print \$1 "." rg ".vcf.gz" }'`
            # Output variant sites
            bcftools concat \${order} -O v > concat.vcf
            vk geno het-polarization concat.vcf > hetpolar.vcf
            bcftools query -f '%CHROM\\t%POS[\\t%GT\\t${SM}\\n]' hetpolar.vcf | grep -v '0/1' | awk -v rg=\${rg} '{ print \$0 "\\t" rg }' > \${rg}.rg_gt.tsv
        done;
        cat *.rg_gt.tsv > rg_gt.tsv
    """
}

process fq_concordance {

    tag{ SM }

    label "R"
    label "xs"

    input:
        tuple file("rg_gt.tsv"), val(SM)

    output:
        file("out_${SM}.tsv")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/fq_concordance.R rg_gt.tsv out_${SM}.tsv
        if [[ ! -e out_${SM}.tsv ]]; then
            touch out_${SM}.tsv
        fi
    """
}

process combine_fq_concordance {

    container null
    executor "local"

    publishDir params.out + "/concordance", mode: 'copy', overwrite: true

    input:
        file("out*.tsv")

    output:
        file("fq_concordance.tsv")

    """
        cat <(echo 'a\tb\tconcordant_sites\ttotal_sites\tconcordance\tSM') out*.tsv > fq_concordance.tsv
    """

}




/* 
    Call variants using the merged site list
*/

process call_variants_union {

    label "vcfkit"
    label "sm"

    tag { SM }

    // stageInMode 'copy'

    input:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi'), path(genome_path), val(genome_basename)

    output:
        val SM, emit: union_vcf_SM
        path "${SM}.union.vcf.gz", emit: union_vcf_set
        path "${SM}.union.vcf.gz.csi", emit: union_vcf_set_indices

    """
        INDEX=`find -L ${genome_path} -name "${genome_basename}.amb" | sed 's/\\.amb\$//'`
        # Re-index the sitelist...because
        # tabix -s 1 -b 2 -e 2 -f sitelist.tsv.gz
        
        contigs="`samtools view -H ${SM}.bam | grep -o -w -e "SN:[A-Z|a-z|0-9|_|-]*" | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs -I {} -P ${task.cpus} sh -c "bcftools mpileup --redo-BAQ -r {} -O b --annotate DP,AD,ADF,ADR,SP --fasta-ref \${INDEX} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels  --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`
        # Output variant sites
        bcftools concat \${order} -O v | \\
        vk geno het-polarization - | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """

}

process generate_union_vcf_list {

    container null
    executor "local"

    publishDir params.out + "/SM", mode: 'copy'

    input:
       val vcf_set

    output:
       file("SM_union_vcfs.txt")

    script:
        print vcf_set

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > SM_union_vcfs.txt
    """
}


process merge_union_vcf_chromosome {

    label "annotation"
    label "sm"

    tag { chrom }

    input:
        tuple file(union_vcfs), val(chrom)

    output:
        tuple val("merge"), file("${chrom}.merged.raw.vcf.gz")

    """
        bcftools merge --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}


process concatenate_union_vcf {

    label "vcfkit"
    label "sm"

    publishDir params.out + "/vcf", mode: 'copy'

    input:
        tuple val("merge"), file("*.merged.raw.vcf.gz"), file("${params.A}.${params.B}.parental.vcf.gz"), file("${params.A}.${params.B}.parental.vcf.gz.tbi")

    output:
        tuple file("NIL.filtered.vcf.gz"), file("NIL.filtered.vcf.gz.csi")

    """
        bcftools concat ${"*.merged.raw.vcf.gz"} | \\
        vk geno het-polarization - | \\
        bcftools view -O z > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
        # Add parental strains
        bcftools merge -O z ${params.A}.${params.B}.parental.vcf.gz merged.raw.vcf.gz > NIL.filtered.vcf.gz
        bcftools index NIL.filtered.vcf.gz
    """
}


process stat_tsv {

    label "annotation"
    label "sm"

    publishDir params.out + "/vcf", mode: 'copy'

    input:
        tuple file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi") 

    output:
        file("NIL.filtered.stats.txt")

    """
        bcftools stats --verbose NIL.filter.vcf.gz > NIL.filtered.stats.txt
    """

}
