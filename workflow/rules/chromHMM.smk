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
        "../envs/sagaconf.yaml"
    threads: config['threads']
    shell:
        'ChromHMM.sh -mx{params.java_max_memory}G BinarizeBam '
        '-gzip -mixed -printposterior '
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
        emissions = "results/chromHMM/learn_model/{num_of_states}/emissions_{num_of_states}.txt",
        overlaps = expand("results/chromHMM/learn_model/{{num_of_states}}/{{num_of_states}}_overlap.txt"),
        posteriors = expand("results/chromHMM/learn_model/{{num_of_states}}/POSTERIOR/{cell_types}_{{num_of_states}}_{chromosome}_posterior.txt", cell_types = cell_types, chromosome = config['chromHMM']['chromosomes'])
    conda:
        "../envs/sagaconf.yaml"
    params:
        chromHMM_jar_loc = config['chromHMM']['chromHMM_jar_loc'],
        java_max_memory = config['chromHMM']['java_max_memory'],
        genome = config['genome']
    threads: config['threads']
    shell:
        'ChromHMM.sh -mx{params.java_max_memory}G LearnModel '
        '-p {threads} '
        'results/chromHMM/binarize_bams '
        'results/chromHMM/learn_model/{wildcards.num_of_states} '
        '{wildcards.num_of_states} '
        '{params.genome}'
