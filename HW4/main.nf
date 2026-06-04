include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'

// 1. Downloading via NCBI
process download_reads {
	tag "$sra_id"
	label 'process_medium'

	input:
		val sra_id
	output:
		tuple val(sra_id), path("${sra_id}_{1,2}.fastq.gz")
	script:
	"""
	fasterq-dump --split-files -e ${task.cpus} ${sra_id}

	gzip ${sra_id}_1.fastq
	gzip ${sra_id}_2.fastq
	"""
}

// 2 & 4. QC
process run_qc {
	tag "${reads_label}_${reads_type}"
	label 'process_low'

	input:
		val reads_type
		tuple val(reads_label), path(reads)
	output:
		path "${reads_label}_${reads_type}_qc_report", type: 'dir'
	script:
	"""

	mkdir ${reads_label}_${reads_type}_qc_report
	fastqc -t ${task.cpus} -o ${reads_label}_${reads_type}_qc_report/ ${reads[0]} ${reads[1]}
	"""
}

// 3. Trimmomatic 
process trimm {
	tag "$reads_label"
	label 'process_medium'

	input:
		tuple val(reads_label), path(reads)
	output:
		tuple val(reads_label), path("${reads_label}_out_R?_p.fq.gz")
	script:
	"""
	export JAVA_TOOL_OPTIONS="-Xmx${task.memory.toGiga()}G"

	
    trimmomatic PE ${reads[0]} ${reads[1]} \\
        ${reads_label}_out_R1_p.fq.gz ${reads_label}_out_R1_u.fq.gz \\
        ${reads_label}_out_R2_p.fq.gz ${reads_label}_out_R2_u.fq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36	
	"""
}

// 5. De novo assembly
process assemble {
    tag "$reads_label"
    label 'process_high'

    input:
        tuple val(reads_label), path(reads)
    output:
        tuple val(reads_label), path("spades_out/scaffolds.fasta"), emit: assembly
    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out --careful -t ${task.cpus}
    """
}

// 6. Mapping reads
process map_reads {
    tag "$reads_label"
    label 'process_high'

    input:
        tuple val(reads_label), path(reference), path(reads)
    output:
        tuple val(reads_label), path("${reads_label}_aligned.sorted.bam"), emit: bam
        tuple val(reads_label), path("${reads_label}_aligned.sorted.bam.bai"), emit: bai
    script:
    """
    bowtie2-build ${reference} ref_index
    
    bowtie2 -p ${task.cpus} -x ref_index -1 ${reads[0]} -2 ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${reads_label}_aligned.sorted.bam
        
    samtools index -@ ${task.cpus} ${reads_label}_aligned.sorted.bam
    """
}

// 7a. Depth calculation
process calculate_depth {
    tag "$reads_label"
    label 'process_low'

    input:
        tuple val(reads_label), path(bam), path(bai)
    output:
        tuple val(reads_label), path("${reads_label}_coverage.tsv"), emit: tsv
    script:
    """
    samtools depth -a ${bam} > ${reads_label}_coverage.tsv
    """
}

// 7b. Plotting
process plot_coverage {
    tag "$reads_label"
    label 'process_low'

    input:
        tuple val(reads_label), path(coverage_tsv)
    output:
        path "${reads_label}_coverage.png", emit: coverage_plot
    script:
    """
    #!/usr/bin/env python3
    import matplotlib.pyplot as plt

    pos, depth = [], []
    with open("${coverage_tsv}") as f:
        for line in f:
            parts = line.split()
            pos.append(int(parts[1]))
            depth.append(int(parts[2]))

    plt.figure(figsize=(12, 4))
    plt.plot(pos, depth, linewidth=0.7)
    plt.xlabel("Position")
    plt.ylabel("Depth")
    plt.title("Coverage for ${reads_label}")
    plt.tight_layout()
    plt.savefig("${reads_label}_coverage.png", dpi=150)
    """
}

process index_fasta {
    tag "$meta.id"
    label 'process_low'
    container 'egoncharov/bowtie2-samtools:2.5.1'

    input:
        tuple val(meta), path(fasta)
    output:
        tuple val(meta), path("${fasta}.fai"), emit: fai
    script:
    """
    samtools faidx ${fasta}
    """
}

process filter_variants {
    tag "$sample_id"
    label 'process_low'
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'

    input:
        tuple val(sample_id), val(group), path(vcf), path(tbi)
    output:
        tuple val(sample_id), val(group), path("${sample_id}.filtered.vcf.gz"), emit: filtered_vcf

    stub:

    """

    echo "STUB: filter_variants for sample ${sample_id} (group: ${group})"
    touch ${sample_id}.filtered.vcf.gz

    """

    script:
    """
    
    bcftools filter \
            --include 'QUAL >= ${params.min_qual} && INFO/DP >= ${params.min_depth}' \
            --output-type z \
            --output ${sample_id}.filtered.vcf.gz \
            ${vcf}
    """

}

workflow trim_and_qc_wf {
    take:
        reads_ch
    main:
        trimm(reads_ch)
        qc_reports = run_qc('trimmed', trimm.out)
    emit:
        trimmed_reads = trimm.out
        qc_reports    = qc_reports
}

workflow reference_wf {
    take:
        trimmed_reads_ch
        reference_path
    main:
        if (reference_path) {
            ref_ch = trimmed_reads_ch.map { label, reads ->
                tuple(label, file(reference_path))
            }
        } else {
            assemble(trimmed_reads_ch)
            ref_ch = assemble.out.assembly
        }
    emit:
        reference = ref_ch
}

workflow mapping_wf {
    take:
        reference_ch
        trimmed_reads_ch
    main:
        mapping_input = reference_ch.join(trimmed_reads_ch)
        map_reads(mapping_input)

        depth_input = map_reads.out.bam.join(map_reads.out.bai)
        calculate_depth(depth_input)

        plot_coverage(calculate_depth.out.tsv)
    emit:
        bam           = map_reads.out.bam
        bai           = map_reads.out.bai
        coverage_plot = plot_coverage.out.coverage_plot
}


// New named workflow for variant calling
workflow variant_calling_wf {
    take:
        bam_ch
        bai_ch
        reference_ch
    main:
        bam_ch
            .join(bai_ch)
            .join(reference_ch)
            .multiMap { label, bam, bai, fasta ->
                bam_input: tuple( [id: label], bam, [], [] )
                ref_input:  tuple( [id: label], fasta )
            }
            .set { split_ch }

        index_fasta(split_ch.ref_input)

        ref_with_fai_ch = split_ch.ref_input
            .join(index_fasta.out.fai)
            .map { meta, fasta, fai -> tuple(meta, fasta, fai) }

        BCFTOOLS_MPILEUP(
            split_ch.bam_input,
            ref_with_fai_ch,
            params.save_mpileup
        )
    emit:
        vcf = BCFTOOLS_MPILEUP.out.vcf
        tbi = BCFTOOLS_MPILEUP.out.tbi
}

workflow {
    main:
        if (params.sample_csv) {
            csv_ch = Channel
                .fromPath(params.sample_csv)
                .splitCsv(header: false)
                .map { row -> tuple(row[0], row[1], file(row[2]), file(row[3])) }

            reads_ch        = csv_ch.map { sample_id, group, r1, r2 -> tuple(sample_id, [r1, r2]) }
            sample_group_ch = csv_ch.map { sample_id, group, r1, r2 -> tuple(sample_id, group) }

            
            reads_ch = csv_ch.map { sample_id, group, r1, r2 -> tuple(sample_id, [r1, r2])}

            sample_group_ch = csv_ch.map { sample_id, group, r1, r2 -> tuple(sample_id, group)}
        }
        else if (params.input_reads_folder) {
            reads_ch = Channel.fromFilePairs(
                "${params.input_reads_folder}/*_{1,2}.{fq,fastq}{,.gz}"
            )

            sample_group_ch = reads_ch.map { label, reads -> tuple(label, 'default')}
        } else {
            download_reads(params.sra_id)
            reads_ch = download_reads.out
            sample_group_ch = reads_ch.map { label, reads -> tuple(label, 'default')}
        }

        initial_qc_res = run_qc('initial', reads_ch)

        trim_and_qc_wf(reads_ch)

        reads_with_group_ch = trim_and_qc_wf.out.trimmed_reads.join(sample_group_ch).map { sample_id, reads, group -> tuple(group, sample_id, reads)}

        split_by_group_ch = reads_with_group_ch.groupTuple(by: 0).flatMap { group, samples, reads_list -> 
                                        [samples, reads_list].transpose().collect { sample, reads -> tuple(sample, reads)}}

        reference_wf(split_by_group_ch, params.reference)

        mapping_wf(reference_wf.out.reference, split_by_group_ch)

        variant_calling_wf(
            mapping_wf.out.bam,
            mapping_wf.out.bai,
            reference_wf.out.reference
        )

        vcf_with_group_ch = variant_calling_wf.out.vcf
            .map { meta, vcf -> tuple(meta.id, vcf) }
            .join(variant_calling_wf.out.tbi.map { meta, tbi -> tuple(meta.id, tbi) })
            .join(sample_group_ch)
            .map { sample_id, vcf, tbi, group -> tuple(sample_id, group, vcf, tbi) }

        filter_variants(vcf_with_group_ch)

        filtered_vcf_ch = filter_variants.out.filtered_vcf

    publish:
        initial_qc      = initial_qc_res
        trimmed_reads   = trim_and_qc_wf.out.trimmed_reads
        trimmed_qc      = trim_and_qc_wf.out.qc_reports
        bam             = mapping_wf.out.bam
        coverage        = mapping_wf.out.coverage_plot
        vcf             = variant_calling_wf.out.vcf
        filtered_vcf    = filtered_vcf_ch
}

output {
    initial_qc      { path 'results/initial_qc' }
    trimmed_reads   { path 'results/trimmed_reads' }
    trimmed_qc      { path 'results/trimmed_qc' }
    bam             { path 'results/bam' }
    coverage        { path 'results/coverage' }
    vcf             { path 'results/vcf' }
    filtered_vcf    { path 'results/filtered_vcf' }
}