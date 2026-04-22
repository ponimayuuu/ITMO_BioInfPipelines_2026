params.input_reads_folder = ''   // path to local fastqs if exists
params.sra_id             = ''   // SRA/DRA accession
params.reference          = ''   // path to a reference, de novo assembly if empty
params.threads            = 20

// 1. Downloading via NCBI
process download_reads {
    conda '/home/egoncharov/miniforge3/envs/tools'
    cpus params.threads

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
    conda '/home/egoncharov/miniforge3/envs/fastqc'
    cpus params.threads

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
    conda '/home/egoncharov/miniforge3/envs/trimmomatic'
    cpus params.threads
    memory '8 GB'

    input:
        tuple val(reads_label), path(reads)
    output:
        tuple val(reads_label), path("${reads_label}_out_R?_p.fq.gz")
    script:
    """
    export JAVA_TOOL_OPTIONS="-Xmx8G"

    trimmomatic PE ${reads[0]} ${reads[1]} \\
        ${reads_label}_out_R1_p.fq.gz ${reads_label}_out_R1_u.fq.gz \\
        ${reads_label}_out_R2_p.fq.gz ${reads_label}_out_R2_u.fq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

// 5. De novo assembly
process assemble {
    conda '/home/egoncharov/miniforge3/envs/spades'
    cpus params.threads

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
    conda '/home/egoncharov/miniforge3/envs/bowtie2'
    cpus params.threads

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
    conda '/home/egoncharov/miniforge3/envs/bowtie2'

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
    conda '/home/egoncharov/miniforge3/envs/python'

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
            ref_ch = trimmed_reads_ch.map { label, reads -> tuple(label, file(reference_path)) }
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
        coverage_plot = plot_coverage.out.coverage_plot
}

workflow {
    main:
        if (params.input_reads_folder) {
            reads_ch = Channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.{fq,fastq}{,.gz}")
        } else {
            download_reads(params.sra_id)
            reads_ch = download_reads.out
        }

        initial_qc_res = run_qc('initial', reads_ch)

        trim_and_qc_wf(reads_ch)

        reference_wf(trim_and_qc_wf.out.trimmed_reads, params.reference)
        mapping_wf(reference_wf.out.reference, trim_and_qc_wf.out.trimmed_reads)

    publish:
        initial_qc    = initial_qc_res
        trimmed_reads = trim_and_qc_wf.out.trimmed_reads
        trimmed_qc    = trim_and_qc_wf.out.qc_reports
        bam           = mapping_wf.out.bam
        coverage      = mapping_wf.out.coverage_plot
}

output {
    initial_qc { path 'results/initial_qc' }
    trimmed_reads { path 'results/trimmed_reads' }
    trimmed_qc { path 'results/trimmed_qc' }
    bam { path 'results/bam' }
    coverage { path 'results/coverage' }
}