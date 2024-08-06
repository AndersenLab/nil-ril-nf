
/*
    ======================
    Generating a site list
    ======================
*/

process generate_sitelist {

    publishDir "${params.out}/sitelist", mode: 'copy'
    
    label "annotation"
    label "sm"

    input:
        file("parental.vcf.gz")

    output:
        tuple file("${params.A}.${params.B}.sitelist.tsv.gz"), file("${params.A}.${params.B}.sitelist.tsv.gz.tbi"), emit: site_list
        tuple val("merge"), file("${params.A}.${params.B}.parental.vcf.gz"), file("${params.A}.${params.B}.parental.vcf.gz.csi"), emit: parental_vcf_only

    """
    # Generate parental VCF
    bcftools view --samples ${params.A},${params.B} -m 2 -M 2 parental.vcf.gz | \\
    awk '\$0 ~ "^#" || (\$0 ~ "0/0" && \$0 ~ "1/1") { print }' | \\
    bcftools filter -O z --include 'FORMAT/GT == "0/0" || FORMAT/GT == "1/1"' > ${params.A}.${params.B}.parental.vcf.gz
    bcftools index ${params.A}.${params.B}.parental.vcf.gz
    # Generate Sitelist
    bcftools query --include 'FORMAT/GT == "0/0" || FORMAT/GT == "1/1"' -f "%CHROM\\t%POS\\t%REF,%ALT\\n" ${params.A}.${params.B}.parental.vcf.gz > ${params.A}.${params.B}.sitelist.tsv
    bgzip ${params.A}.${params.B}.sitelist.tsv
    tabix -s 1 -b 2 -e 2 ${params.A}.${params.B}.sitelist.tsv.gz
    """   
}


/*
    ===============
    Fastq alignment
    ===============
*/

process perform_alignment {

    tag { ID }

    label "trim"
    label "sm"

    input:
        tuple val(SM), val(ID), val(LB), file(fq1), file(fq2), path(genome_path), val(genome_basename)

    output:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai"), emit: aligned_bams
        tuple val(SM), file("${ID}.bam"), file("${ID}.bam.bai"), emit: sample_aligned_bams

    """
        INDEX=`find -L ${genome_path} -name "${genome_basename}.amb" | sed 's/\\.amb\$//'`
        bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' \${INDEX} ${fq1} ${fq2} | \\
        samtools view -hb -@ ${task.cpus} > mapped.bam
        samtools sort -@ ${task.cpus} -T ${params.tmpdir}/${ID} -o ${ID}.bam mapped.bam
        rm mapped.bam
        samtools index -@ ${task.cpus} ${ID}.bam
        if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}

/*
    fq idx stats
*/

process fq_idx_stats {
    
    tag { ID }

    label "trim"
    label "sm"

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")

    output:
        file 'fq_idxstats'

    """
        samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {
    
    container null
    executor "local"

    publishDir params.out + "/fq", mode: 'copy'

    input:
        val bam_idxstats

    output:
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> fq_bam_idxstats.tsv
    """

}
