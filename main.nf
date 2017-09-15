#!/usr/bin/env nextflow
/*
    Set these parameters in nextflow.config
*/

/* 
 * Authors: 
 * - Daniel Cook <danielecook@gmail.com>
 */

date = new Date().format( 'yyyy-MM-dd' )

/*
    =============
    Configuration
    =============
*/
debug="${params.debug?:'false'}"

// Define parents
A="WN2001"
B="${params.B?:'CB4856'}"

// Analysis Directory
analysis_dir = params.analysis_prefix + "/NIL-${A}-${B}-${params.analysis_name?:date}"

// Reference
File f = new File(params.reference);
reference = f.getAbsolutePath();

// Define VCF
parental_vcf=file("${params.vcf}")
hmm_plot_script=Channel.fromPath("${workflow.projectDir}/plot_hmm.R")
cross_object_script=file("${workflow.projectDir}/generate_cross_object.R")
qc_plots_script=file("${workflow.projectDir}/qc_plots.R")

if (!parental_vcf.exists()) {
    println("You must define a parental vcf with --vcf")
    System.exit(0)
}

if (params.debug == 'true') {
    println "Using debug mode"
}

// Define contigs here!
contig_list = ["I", "II", "III", "IV", "V", "X", "MtDNA"];
contigs = Channel.from(contig_list)

println "Processing NIL Sequence data located in ${params.fq_folder}"
println "Using Reference: ${reference}" 
println "Using Parents GT: ${A} - ${B}"
println "Outputting results to ${analysis_dir}"

strainFile = new File(params.fq_folder + "/fq_sheet.tsv")
fqs = Channel.from(strainFile.collect { it.tokenize( '\t' ) })

/*
    ======================
    Generating a site list
    ======================
*/

process generate_sitelist {


    publishDir analysis_dir + "/sitelist", mode: 'copy'
    
    input:
        file("parental.vcf.gz") from parental_vcf

    output:
        set file("${A}.${B}.sitelist.tsv.gz"), file("${A}.${B}.sitelist.tsv.gz.tbi") into site_list
        set file("${A}.${B}.parental.vcf.gz"), file("${A}.${B}.parental.vcf.gz.csi") into parental_vcf_only

    """
        # Generate parental VCF
        bcftools view --samples ${A},${B} -m 2 -M 2 ${parental_vcf} | \
        vk filter ALT --min=1 - | \
        vk filter REF --min=1 - | \
        bcftools filter -O z --include 'FORMAT/GT == "0/0" || FORMAT/GT == "1/1"' > ${A}.${B}.parental.vcf.gz
        bcftools index ${A}.${B}.parental.vcf.gz

        # Generate Sitelist
        bcftools query --include 'FORMAT/GT == "0/0" || FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' ${A}.${B}.parental.vcf.gz| \
        bgzip -c > ${A}.${B}.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 ${A}.${B}.sitelist.tsv.gz
    """   
}


/*
    ===============
    Fastq alignment
    ===============
*/

process perform_alignment {

    echo true

    cpus params.align_cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2 from fqs
    output:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") into aligned_bams
        set SM, file("${ID}.bam") into sample_aligned_bams

    script:
    if(debug == 'true')
        """
            zcat ${fq1} | head -n 5000 | gzip > fq1.fq.gz
            zcat ${fq2} | head -n 5000 | gzip > fq2.fq.gz
            bwa mem -t ${params.align_cores} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference} fq1.fq.gz fq2.fq.gz | \\
            sambamba view --nthreads=${params.align_cores} --sam-input --format=bam --with-header /dev/stdin | \\
            sambamba sort --nthreads=${params.align_cores} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
            sambamba index --nthreads=${params.align_cores} ${ID}.bam

            if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
                exit 1;
            fi
        """ 
    else
        """
            bwa mem -t ${params.align_cores} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference} ${fq1} ${fq2} | \\
            sambamba view --nthreads=${params.align_cores} --sam-input --format=bam --with-header /dev/stdin | \\
            sambamba sort --nthreads=${params.align_cores} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
            sambamba index --nthreads=${params.align_cores} ${ID}.bam

            if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
                exit 1;
            fi
        """
}

aligned_bams.into { 
                           sample_bams_fq_idx_stats;
                           fq_stat_bams;
                           fq_cov_bam_indices; 
                         }

/*
    fq idx stats
*/

process fq_idx_stats {
    
    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from sample_bams_fq_idx_stats
    output:
        file 'fq_idxstats' into fq_idxstats_set

    """
        samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val bam_idxstats from fq_idxstats_set.toSortedList()

    output:
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> fq_bam_idxstats.tsv
    """

}

/*
    fq bam stats
*/

process fq_bam_stats {

    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from fq_stat_bams

    output:
        file 'bam_stat' into bam_stat_files

    """
        cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_bam_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val stat_files from bam_stat_files.toSortedList()

    output:
        file("fq_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
        cat ${stat_files.join(" ")} >> fq_bam_stats.tsv
    """
}

/*
    Fastq coverage
*/
process fq_coverage {

    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from fq_cov_bam_indices
    output:
        file("${ID}.coverage.tsv") into fq_coverage


    """
        bam coverage ${ID}.bam > ${ID}.coverage.tsv
    """
}

process fq_coverage_merge {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val fq_set from fq_coverage.toSortedList()

    output:
        file("fq_coverage.full.tsv")
        file("fq_coverage.tsv") into fq_coverage_plot

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}

sample_aligned_bams.into { sample_aligned_bams_out; sample_aligned_bams_use}

sample_aligned_bams_out.groupTuple().println()

process merge_bam {

    echo true

    storeDir analysis_dir + "/bam"

    cpus params.align_cores

    tag { SM }

    input:
        set SM, bam from sample_aligned_bams_use.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into merged_SM
        file("${SM}.duplicates.txt") into duplicates_file

    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam[0]} ${SM}.merged.bam
        ln -s ${bam[0]}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${params.align_cores} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${params.align_cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${params.align_cores} ${SM}.bam
    """
}

merged_SM.into { 
    bams_stats;
    bams_idxstats;
    merged_bams_for_coverage;
    merged_bams_union;
}



/*
 SM idx stats
*/

process idx_stats_SM {
    
    tag { SM }

    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
    output:
        file 'SM_bam_idxstats' into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > SM_bam_idxstats
    """
}

process combine_idx_stats {

    publishDir analysis_dir +"/SM", mode: 'copy'

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv") into SM_bam_idxstats_plot

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
    """

}


/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    publishDir analysis_dir + "/SM", mode: 'copy'

    input:
        val stat_files from SM_bam_stat_files.toSortedList()

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}



process format_duplicates {

    publishDir analysis_dir + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv") into bam_duplicates_plot

    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}


/*
    Coverage Bam
*/
process SM_coverage {

    tag { SM }

    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into SM_coverage

    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}


process SM_coverage_merge {

    publishDir analysis_dir + "/SM", mode: 'copy'


    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv") into SM_coverage_plot
        file("SM_coverage.tsv") into SM_coverage_cross_obj

    """
        echo -e 'SM\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
    """

}


/* 
    Call variants using the merged site list
*/

merged_bams_union.merge(site_list) { a, b -> [a, b].flatten() }.into { union_vcf_channel; union_vcf_channel_print }

union_vcf_channel_print.println()



/*
    Call variants
*/

process call_variants_union {

    echo true

    cpus params.variant_cores

    tag { SM }

    stageInMode 'copy'

    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_channel

    output:
        val SM into union_vcf_SM
        file("${SM}.union.vcf.gz") into union_vcf_set
        file("${SM}.union.vcf.gz.csi") into union_vcf_set_indices

    """
        # Re-index the sitelist...because
        # tabix -s 1 -b 2 -e 2 -f sitelist.tsv.gz
        
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${params.variant_cores} sh -c "samtools mpileup --redo-BAQ --region {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${params.reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels  --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
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

    cpus 1 

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
       val vcf_set from union_vcf_set.toSortedList()

    output:
       file("union_vcfs.txt") into union_vcfs

    script:
        print vcf_set

    """
        # run 3
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
    """
}

union_vcfs_in = union_vcfs.spread(contigs)

process merge_union_vcf_chromosome {

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in

    output:
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        bcftools merge --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}

// Generate a list of ordered files.
contig_raw_vcf = contig_list*.concat(".merged.raw.vcf.gz")

process concatenate_union_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        val merge_vcf from raw_vcf.toSortedList()
        set file("${A}.${B}.parental.vcf.gz"), file("${A}.${B}.parental.vcf.gz.tbi") from parental_vcf_only

    output:
        set file("NIL.filtered.vcf.gz"), file("NIL.filtered.vcf.gz.csi") into filtered_vcf

    """
        # First - concatenate chrom sets.
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat ${contig_raw_vcf.join(" ")} | \\
        vk geno het-polarization - | \\
        bcftools view -O z > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz

        # Add parental strains
        bcftools merge -O z ${A}.${B}.parental.vcf.gz merged.raw.vcf.gz > NIL.filtered.vcf.gz
        bcftools index NIL.filtered.vcf.gz
    """
}


filtered_vcf.into { filtered_vcf_stat; hmm_vcf; hmm_vcf_clean; hmm_vcf_out; vcf_tree }


process stat_tsv {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi")  from filtered_vcf_stat

    output:
        file("NIL.filtered.stats.txt")

    """
        bcftools stats --verbose NIL.filter.vcf.gz > NIL.filtered.stats.txt
    """

}

process output_hmm {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi") from hmm_vcf

    output:
        file("gt_hmm.tsv")

    """
        {
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        } || { 
            echo 'is yahmm installed?' 
        }
        vk hmm --A=${A} --B=${B} NIL.filter.vcf.gz > gt_hmm.tsv
    """

}

process output_hmm_clean {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("NIL.filter.vcf.gz"), file("NIL.filter.vcf.gz.csi") from hmm_vcf_clean

    output:
        file("gt_hmm_fill.tsv") into gt_hmm_fill

    """
        {
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        } || { 
            echo 'is yahmm installed?' 
        }
        vk hmm --transition=1e-12 --infill --endfill --A=${A} --B=${B} NIL.filter.vcf.gz > gt_hmm_fill.tsv
    """

}


process output_hmm_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("NIL.vcf.gz"), file("NIL.vcf.gz.csi") from hmm_vcf_out

    output:
        set file("NIL.hmm.vcf.gz"), file("NIL.hmm.vcf.gz.csi") into gt_hmm

    """
        {
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        } || { 
            echo 'is yahmm installed?' 
        }
        vk hmm --transition=1e-12 --vcf-out --all-sites --A=${A} --B=${B} NIL.vcf.gz | bcftools view -O z > NIL.hmm.vcf.gz
        bcftools index NIL.hmm.vcf.gz
    """

}

gt_hmm.into { gt_hmm_tsv; gt_hmm_vcf }

process plot_hmm {

    publishDir analysis_dir + "/hmm", mode: 'copy'

    input:
        file("gt_hmm_fill.tsv") from gt_hmm_fill
        file("script.R") from hmm_plot_script

    output:
        file("gt_hmm.png")
        file("gt_hmm.svg")

    """
        Rscript --vanilla script.R
    """

}

process generate_issue_plots {

    publishDir analysis_dir + "/plots", mode: 'copy'

    errorStrategy 'ignore'

    input:
        file("SM_bam_idxstats.tsv") from SM_bam_idxstats_plot
        file("fq_coverage.tsv") from fq_coverage_plot
        file("SM_coverage.tsv") from SM_coverage_plot
        file("bam_duplicates.tsv") from bam_duplicates_plot
        file("qc_plots.R") from qc_plots_script

    output:
        file("coverage_comparison.png")
        file("coverage_comparison.svg")
        file("unmapped_reads.png")
        file("unmapped_reads.svg")
        file("duplicates.png")
        file("duplicates.svg")

    """
        Rscript --vanilla qc_plots.R
    """
}


process output_tsv {

    publishDir analysis_dir + "/hmm", mode: 'copy'

    input:
        set file("NIL.hmm.vcf.gz"), file("NIL.hmm.vcf.gz.csi") from gt_hmm_tsv

    output:
        file("gt_hmm.tsv")

    """
        cat <(echo -e "CHROM\tPOS\tSAMPLE\tGT\tGT_ORIG\tAD") <(bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\t%GT_ORIG\t%AD\n]' NIL.hmm.vcf.gz | sed 's/0\\/0/0/g' | sed 's/1\\/1/1/g') > gt_hmm.tsv
    """

}
