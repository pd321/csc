rule binarize_bams:
    input:
        metadata_file = "config/chromHMM_metadata.tsv",
        genome_sizes = config['genome_sizes']
    output:
        binarized_bams = expand("results/chromHMM/binarize_bams/{cell_types}_{chromosome}_binary.txt.gz",
            chromosome = config['chromHMM']['chromosomes'], cell_types = cell_types)
    params:
        bam_directory = config['chromHMM']['bam_directory'],
        java_max_memory = config['chromHMM']['java_max_memory']
    conda:
        "../envs/chromHMM.yaml"
    threads: config['threads']
    shell:
        'ChromHMM.sh -Xmx{params.java_max_memory}G BinarizeBam '
        '-gzip -mixed '
        '{input.genome_sizes} '
        '{params.bam_directory} '
        '{input.metadata_file} '
        'results/chromHMM/binarize_bams'

# Learn the model from binarized bams
rule learn_model:
    input:
        binarized_bams = rules.binarize_bams.output.binarized_bams
    output:
        model = "results/chromHMM/learn_model/{num_of_states}/model_{num_of_states}.txt",
        emissions = "results/chromHMM/learn_model/{num_of_states}/emissions_{num_of_states}.txt"
    conda:
        "../envs/chromHMM.yaml"
    params:
        java_max_memory = config['chromHMM']['java_max_memory'],
        genome = config['genome']
    threads: config['threads']
    shell:
        'ChromHMM.sh -Xmx{params.java_max_memory}G LearnModel '
        '-p {threads} '
        'results/chromHMM/binarize_bams '
        'results/chromHMM/learn_model/{wildcards.num_of_states} '
        '{wildcards.num_of_states} '
        '{params.genome}'
