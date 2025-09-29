/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNGSPipeline.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.samples, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.samples) { ch_samples = file(params.samples) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_dummy_file              = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { QC_FASTQ          } from '../subworkflows/local/qc_fastq'
include { ALIGN_STAR        } from '../subworkflows/local/align_star'
include { ALIGN_BOWTIE2     } from '../subworkflows/local/align_bowtie2'
include { ALIGN_BISMARK     } from '../subworkflows/local/align_bismark'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { QUANTIFY_GENES    } from '../subworkflows/local/quantify_genes'
include { DIFFERENTIAL_DESEQ2 } from '../subworkflows/local/differential_deseq2'
include { CALL_PEAKS        } from '../subworkflows/local/call_peaks'
include { CALL_VARIANTS     } from '../subworkflows/local/call_variants'
include { METHYLATION_ANALYSIS } from '../subworkflows/local/methylation_analysis'
include { ML_INTEGRATION    } from '../subworkflows/local/ml_integration'
include { PATHWAY_ANALYSIS  } from '../subworkflows/local/pathway_analysis'
include { VISUALIZATION     } from '../subworkflows/local/visualization'
include { FINAL_REPORTING   } from '../subworkflows/local/final_reporting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ngs_pipeline {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_samples
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Prepare genome references and indices
    //
    PREPARE_GENOME (
        params.genome_fasta,
        params.genome_gtf,
        params.star_index,
        params.bowtie2_index,
        params.bismark_index
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // Branch samples by assay type for appropriate processing
    //
    INPUT_CHECK.out.reads
        .branch { meta, reads ->
            rnaseq: meta.assay == 'RNAseq' || meta.assay == 'EUseq'
            sdripseq: meta.assay == 'sDRIPseq' 
            endseq: meta.assay == 'ENDseq'
            bsseq: meta.assay == 'BSseq'
        }
        .set { ch_reads_branched }

    //
    // SUBWORKFLOW: FastQC and adapter trimming
    //
    QC_FASTQ (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(QC_FASTQ.out.versions)

    //
    // SUBWORKFLOW: Align RNA-seq and EU-seq with STAR
    //
    ALIGN_STAR (
        ch_reads_branched.rnaseq,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.gtf
    )
    ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

    //
    // SUBWORKFLOW: Align sDRIP-seq and END-seq with Bowtie2
    //
    ch_bowtie2_input = ch_reads_branched.sdripseq.mix(ch_reads_branched.endseq)
    
    ALIGN_BOWTIE2 (
        ch_bowtie2_input,
        PREPARE_GENOME.out.bowtie2_index
    )
    ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions)

    //
    // SUBWORKFLOW: Align Bisulfite-seq with Bismark
    //
    ALIGN_BISMARK (
        ch_reads_branched.bsseq,
        PREPARE_GENOME.out.bismark_index
    )
    ch_versions = ch_versions.mix(ALIGN_BISMARK.out.versions)

    //
    // Combine all aligned BAM files
    //
    ch_aligned_bams = ALIGN_STAR.out.bam
        .mix(ALIGN_BOWTIE2.out.bam)
        .mix(ALIGN_BISMARK.out.bam)

    //
    // SUBWORKFLOW: Sort BAMs and collect alignment statistics
    //
    BAM_SORT_STATS_SAMTOOLS (
        ch_aligned_bams,
        PREPARE_GENOME.out.fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // SUBWORKFLOW: Gene quantification for RNA-seq/EU-seq
    //
    ch_rnaseq_bams = BAM_SORT_STATS_SAMTOOLS.out.bam
        .filter { meta, bam, bai -> meta.assay == 'RNAseq' || meta.assay == 'EUseq' }

    QUANTIFY_GENES (
        ch_rnaseq_bams,
        PREPARE_GENOME.out.gtf
    )
    ch_versions = ch_versions.mix(QUANTIFY_GENES.out.versions)

    //
    // SUBWORKFLOW: Differential expression analysis
    //
    if (!params.skip_deseq2) {
        DIFFERENTIAL_DESEQ2 (
            QUANTIFY_GENES.out.counts,
            INPUT_CHECK.out.samplesheet
        )
        ch_versions = ch_versions.mix(DIFFERENTIAL_DESEQ2.out.versions)
    }

    //
    // SUBWORKFLOW: Peak calling for sDRIP-seq and END-seq
    //
    ch_peak_calling_bams = BAM_SORT_STATS_SAMTOOLS.out.bam
        .filter { meta, bam, bai -> meta.assay == 'sDRIPseq' || meta.assay == 'ENDseq' }

    CALL_PEAKS (
        ch_peak_calling_bams
    )
    ch_versions = ch_versions.mix(CALL_PEAKS.out.versions)

    //
    // SUBWORKFLOW: Variant calling (optional)
    //
    if (params.enable_variants) {
        CALL_VARIANTS (
            ch_rnaseq_bams,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fasta_fai,
            PREPARE_GENOME.out.dict,
            params.dbsnp_vcf
        )
        ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)
    }

    //
    // SUBWORKFLOW: Methylation analysis for Bisulfite-seq
    //
    ch_bismark_bams = BAM_SORT_STATS_SAMTOOLS.out.bam
        .filter { meta, bam, bai -> meta.assay == 'BSseq' }

    if (ch_bismark_bams) {
        METHYLATION_ANALYSIS (
            ch_bismark_bams,
            PREPARE_GENOME.out.fasta
        )
        ch_versions = ch_versions.mix(METHYLATION_ANALYSIS.out.versions)
    }

    //
    // SUBWORKFLOW: Visualization tracks and sessions
    //
    VISUALIZATION (
        BAM_SORT_STATS_SAMTOOLS.out.bam,
        CALL_PEAKS.out.peaks.collect().ifEmpty([]),
        PREPARE_GENOME.out.gtf
    )
    ch_versions = ch_versions.mix(VISUALIZATION.out.versions)

    //
    // SUBWORKFLOW: Machine learning integration
    //
    if (!params.skip_ml) {
        ML_INTEGRATION (
            QUANTIFY_GENES.out.counts.collect().ifEmpty([]),
            CALL_PEAKS.out.peaks.collect().ifEmpty([]),
            INPUT_CHECK.out.samplesheet
        )
        ch_versions = ch_versions.mix(ML_INTEGRATION.out.versions)
    }

    //
    // SUBWORKFLOW: Pathway and network analysis
    //
    if (params.enable_string || params.enable_tcga) {
        PATHWAY_ANALYSIS (
            DIFFERENTIAL_DESEQ2.out.results.collect().ifEmpty([]),
            params.enable_string,
            params.enable_tcga
        )
        ch_versions = ch_versions.mix(PATHWAY_ANALYSIS.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNGSPipeline.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNGSPipeline.methodsDescriptionText(workflow, ch_multiqc_custom_config, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(QC_FASTQ.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_FASTQ.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_BOWTIE2.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_BISMARK.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_GENES.out.summary.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    //
    // SUBWORKFLOW: Final consolidated reporting
    //
    FINAL_REPORTING (
        multiqc_report,
        DIFFERENTIAL_DESEQ2.out.results.collect().ifEmpty([]),
        CALL_PEAKS.out.peaks.collect().ifEmpty([]),
        ML_INTEGRATION.out.results.collect().ifEmpty([]),
        PATHWAY_ANALYSIS.out.results.collect().ifEmpty([]),
        VISUALIZATION.out.sessions.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix(FINAL_REPORTING.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}