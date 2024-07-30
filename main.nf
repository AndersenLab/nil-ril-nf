#!/usr/bin/env nextflow
/* 
 * Authors: 
 * - Daniel Cook <danielecook@gmail.com>
 * - Katie Evans <katiesevans9@gmail.com>
 *  
 */

 if( !nextflow.version.matches('23.0+') ) {
    println "This workflow requires Nextflow version 23.0 or greater -- You are running version $nextflow.version"
    println "On Rockfish, you can use `module load python/anaconda; source activate /data/eande106/software/conda_envs/nf23_env`"
    //exit 1
}

nextflow.enable.dsl=2


params.out = "NIL-${params.A}-${params.B}-${params.day}"
date = new Date().format( 'yyyyMMdd' )

// debug
if (params.debug == true) {
    println """
        ***Using debug mode***
    """
    params.relative = true
    params.fqs = "${workflow.projectDir}/test_data/fq_sheet.tsv"
    params.vcf = "${workflow.projectDir}/test_data/N2_CB.simple.vcf.gz"
    params.reference = "${params.data_dir}/c_elegans/genomes/PRJNA13758/${params.genome}/c_elegans.PRJNA13758.${params.genome}.genome.fa.gz"
} else {
    params.relative = false
    params.fqs = "(required)"
    params.vcf = "(required)"
    params.reference = "(required)"
}

// checks
if (params.vcf == "(required)" || params.reference == "(required)" || params.fqs == "(required)") {
    println """

    Error: VCF, Reference, and FQ sheet are required for analysis. Please check parameters.

    VCF: ${params.vcf}
    Reference: ${params.reference}
    FQs: ${params.fqs}

    """
    System.exit(1)
} 


param_summary = '''

███╗   ██╗██╗██╗      ██████╗ ██╗██╗      ███╗   ██╗███████╗
████╗  ██║██║██║      ██╔══██╗██║██║      ████╗  ██║██╔════╝
██╔██╗ ██║██║██║█████╗██████╔╝██║██║█████╗██╔██╗ ██║█████╗  
██║╚██╗██║██║██║╚════╝██╔══██╗██║██║╚════╝██║╚██╗██║██╔══╝  
██║ ╚████║██║███████╗ ██║  ██║██║███████╗ ██║ ╚████║██║     
╚═╝  ╚═══╝╚═╝╚══════╝ ╚═╝  ╚═╝╚═╝╚══════╝ ╚═╝  ╚═══╝╚═╝     
                                                            
''' + """
    parameters           description                    Set/Default
    ==========           ===========                    =======
    
    --debug              Set to 'true' to test          ${params.debug}
    --A                  Parent A                       ${params.A}
    --B                  Parent B                       ${params.B}
    --cA                 Parent A color (for plots)     ${params.cA}
    --cB                 Parent B color (for plots)     ${params.cB}
    --out                Directory to output results    ${params.out}
    --fqs                fastq file (see help)          ${params.fqs}
    --relative           use relative fastq prefix      ${params.relative}
    --reference          Reference Genome               ${params.reference}
    --vcf                VCF to fetch parents from      ${params.vcf}
    --transition         Transition Prob                ${params.transition}
    --tmpdir             A temporary directory          ${params.tmpdir}
    --email              Email to be sent results       ${params.email}
    HELP: http://andersenlab.org/dry-guide/pipeline-nil/
  ---------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
"""

println param_summary
if (params.help) {
    exit 1
}

// Includes
include {kmer_counting; merge_kmer} from './modules/kmers.nf'
include {generate_sitelist; perform_alignment; fq_idx_stats; fq_combine_idx_stats} from './modules/alignment.nf'
include {fq_bam_stats; combine_bam_stats; combine_SM_bam_stats} from './modules/stats.nf'
include {idx_stats_SM; combine_idx_stats; SM_bam_stats} from './modules/stats.nf'
include {merge_bam; format_duplicates} from './modules/bam.nf'
include {fq_coverage; fq_coverage_merge; SM_coverage; SM_coverage_merge} from './modules/coverage.nf'
include {split_fq; fq_concordance; combine_fq_concordance} from './modules/variants.nf'
include {call_variants_union; generate_union_vcf_list; stat_tsv} from './modules/variants.nf'
include {concatenate_union_vcf; merge_union_vcf_chromosome} from './modules/variants.nf'
include {output_hmm; output_hmm_fill; output_hmm_vcf; plot_hmm} from './modules/hmm.nf'
include {generate_issue_plots; output_tsv; generate_cross_object} from './modules/outputs.nf'

// Generate workflow
workflow {

    genome_path = "$params.reference".substring(0, "$params.reference".lastIndexOf("/"))
    genome_basename = "$params.reference".substring("$params.reference".lastIndexOf("/") + 1)

    if (params.relative) {
        fqs = Channel.fromPath(params.fqs, checkIfExists: true)
            .ifEmpty { exit 1, "sample sheet not found" }
            .splitCsv(sep: '\t')
            .map { SM, ID, LB, fq1, fq2 -> [
                SM, ID, LB, file("${workflow.projectDir}/${fq1}"),
                file("${workflow.projectDir}/${fq2}"), genome_path,
                genome_basename] }
    } else {
        fqs = Channel.fromPath(params.fqs, checkIfExists: true)
            .ifEmpty { exit 1, "sample sheet not found" }
            .splitCsv(sep: '\t')
            .map { SM, ID, LB, fq1, fq2 -> [
                SM, ID, LB, file(fq1), file(fq2), genome_path, genome_basename] }
    }

    Channel.fromFile(params.vcf, checkIfExists: true)
    Channel.fromFile(params.reference, checkIfExists: true)

    // generate site list
    Channel.fromPath(params.vcf) | generate_sitelist

    // kmers
    fqs | kmer_counting
    kmer_counting.out.collect() | merge_kmer

    // alignment
    fqs | perform_alignment
    perform_alignment.out.sample_aligned_bams.groupTuple() | merge_bam
    merge_bam.out.merged_SM | SM_coverage
    SM_coverage.out.toSortedList() | SM_coverage_merge
    perform_alignment.out.aligned_bams
        .combine(generate_sitelist.out.site_list)
        .combine(Channel.fromPath(genome_path))
        .combine(Channel.of(genome_basename)) | split_fq
    split_fq.out | fq_concordance
    fq_concordance.out
        .toSortedList() | combine_fq_concordance
    perform_alignment.out.aligned_bams
    split_fq.out | fq_coverage
    fq_coverage.out
        .toSortedList() | fq_coverage_merge
    merge_bam.out.duplicates_file
        .toSortedList() | format_duplicates

    // call variants
    // merge_bam.out.merged_SM.combine(generate_sitelist.out.site_list) | call_variants_union
    // call_variants_union.out.union_vcf_set.toSortedList() | generate_union_vcf_list
    // generate_union_vcf_list.out
    //     .spread(Channel.of(["I", "II", "III", "IV", "V", "X", "MtDNA"])) | merge_union_vcf_chromosome
    // merge_union_vcf_chromosome.out
    //     .groupTuple() 
    //     .join(generate_sitelist.out.parental_vcf_only) | concatenate_union_vcf

    // // stats
    // merge_bam.out.merged_SM | SM_bam_stats
    // SM_bam_stats.out
    //     .toSortedList() | combine_SM_bam_stats
    // merge_bam.out.merged_SM | idx_stats_SM
    // idx_stats_SM.out
    //     .toSortedList() | combine_idx_stats
    // concatenate_union_vcf.out | stat_tsv
    // perform_alignment.out.aligned_bams | fq_idx_stats
    // fq_idx_stats.out
    //     .toSortedList() | fq_combine_idx_stats
    // perform_alignment.out.aligned_bams | fq_bam_stats
    // fq_bam_stats.out
    //     .toSortedList() | combine_bam_stats

    // // hmm
    // concatenate_union_vcf.out | output_hmm
    // concatenate_union_vcf.out | output_hmm_fill | plot_hmm
    // concatenate_union_vcf.out | output_hmm_vcf | output_tsv

    // // plot issues
    // fq_coverage_merge.out.fq_coverage_plot
    //     .combine(combine_idx_stats.out)
    //     .combine(SM_coverage_merge.out.SM_coverage_plot)
    //     .combine(format_duplicates.out) | generate_issue_plots 

    // // generate cross object geno for RILs
    // if(params.cross_obj) {
    //     output_hmm_vcf.out
    //         .combine(SM_coverage_merge.out.SM_coverage_plot)
    //         .combine(output_hmm_fill.out) | generate_cross_object
    // }


}

workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    Parameters
    ----------
    debug              Set to 'true' to test          ${params.debug}
    A                  Parent A                       ${params.A}
    B                  Parent B                       ${params.B}
    cA                 Parent A color (for plots)     ${params.cA}
    cB                 Parent B color (for plots)     ${params.cB}
    out                Directory to output results    ${params.out}
    fqs                fastq file (see help)          ${params.fqs}
    relative           use relative fastq prefix      ${params.relative}
    reference          Reference Genome               ${reference_handle}
    vcf                VCF to fetch parents from      ${params.vcf}
    transition         Transition Prob                ${params.transition}
    """

    println summary

    // mail summary
    //['mail', '-s', 'nil-ril-nf', params.email].execute() << summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }



}