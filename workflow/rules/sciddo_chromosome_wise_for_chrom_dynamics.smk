# This is needed for the extraction of specific changes like regions that go from e.g. state 9 or 8 to state 1 or 2
# If we don't need the heatmap, we can skip this step
rule convert_states_and_score_by_chromosome:
    input:
        model_emissions="results/chromHMM/learn_model/{sciddo_chromHMM_model}/model_{sciddo_chromHMM_model}.txt",
        chrom_sizes=config["genome_sizes"],
        state_labels="config/sciddo_state_labels_for_state_{sciddo_chromHMM_model}.tsv",
        state_colours="config/sciddo_state_colours_for_state_{sciddo_chromHMM_model}.tsv",
        design_matrix="config/sciddo_design_matrix_for_state_{sciddo_chromHMM_model}.tsv",
    output:
        sciddo_convert_out="results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/convert/sciddo_convert-{chromosome}.h5",
    conda:
        "../envs/sciddo.yaml"
    threads: config["threads"]
    shell:
        "sciddo.py "
        "--no-conf-dump "
        "--workers {threads} "
        "convert "
        "--state-seg results/chromHMM/learn_model/{wildcards.sciddo_chromHMM_model}/ "
        "--chrom-sizes {input.chrom_sizes} "
        "--state-labels {input.state_labels} "
        "--state-colors {input.state_colours} "
        "--design-matrix {input.design_matrix} "
        "--model-emissions {input.model_emissions} "
        "--output {output.sciddo_convert_out} "
        '--chrom-filter "{wildcards.chromosome}$" '
        "&& "
        "sciddo.py "
        "--no-conf-dump "
        "--workers {threads} "
        "stats "
        "--counts "
        "--agreement "
        "--sciddo-data {output.sciddo_convert_out} "
        "&& "
        "sciddo.py "
        "--no-conf-dump "
        "score "
        "--sciddo-data {output.sciddo_convert_out} "
        "--add-scoring emission "
        "--treat-background penalized"

rule scan_by_chromosome:
    input:
        sciddo_convert_out=rules.convert_states_and_score_by_chromosome.output.sciddo_convert_out,
    output:
        sciddo_scan_out="results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/scan/{sciddo_group1}-vs-{sciddo_group2}-{chromosome}.h5",
    conda:
        "../envs/sciddo.yaml"
    threads: config["threads"]
    shell:
        "sciddo.py "
        "--no-conf-dump "
        "--workers {threads} "
        "scan "
        "--sciddo-data {input.sciddo_convert_out} "
        "--run-out {output.sciddo_scan_out} "
        "--scoring penem "
        "--adjust-group-length linear "
        "--merge-segments "
        "--compute-raw-stats 0 "
        "--compute-merged-stats 0 "
        "--select-groups "
        "--group1 {wildcards.sciddo_group1} "
        "--group2 {wildcards.sciddo_group2}"

rule dump_chromatin_dynamics:
    input:
        sciddo_convert_out=rules.convert_states_and_score_by_chromosome.output.sciddo_convert_out,
        sciddo_scan_out=rules.scan_by_chromosome.output.sciddo_scan_out,
    output:
        sciddo_dump_out="results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/dump/{from_state}-to-{to_state}/{sciddo_group1}-vs-{sciddo_group2}-{chromosome}.tsv",
    conda:
        "../envs/sciddo.yaml"
    threads: config["threads"]
    params:
        from_state=lambda wildcards: config["sciddo"]["chromatin_dynamics"][
            wildcards.from_state
        ],
        to_state=lambda wildcards: config["sciddo"]["chromatin_dynamics"][
            wildcards.to_state
        ],
    shell:
        "sciddo.py "
        "--no-conf-dump "
        "dump "
        "--data-file {input.sciddo_convert_out} "
        "--support-file {input.sciddo_scan_out} "
        "--output {output.sciddo_dump_out} "
        "--data-type dynamics "
        "--from-states {params.from_state} "
        "--to-states {params.to_state} "
        "--split-segments --keep-duplicates "
        " || touch {output.sciddo_dump_out}"

rule dump_chromatin_dynamics_to_bed:
    input:
        sciddo_dump_out=rules.dump_chromatin_dynamics.output.sciddo_dump_out,
    output:
        sciddo_dump_bed_out="results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/dump-bed/{from_state}-to-{to_state}/{sciddo_group1}-vs-{sciddo_group2}-{chromosome}.bed",
    params:
        req_in_how_many_replicates=config["sciddo"]["chromatin_dynamics"][
            "req_in_how_many_replicates"
        ],
    conda:
        "../envs/bedtools.yaml"
    shell:
        "tail -n +2 {input.sciddo_dump_out} | "
        "cut -f1,2,3 | "
        "sort | uniq -c | "
        "awk ' {{ if ( $1+0 >= {params.req_in_how_many_replicates} ) print $2,$3,$4 }} ' OFS='\\t' | "
        "bedtools sort > {output.sciddo_dump_bed_out}"


rule merge_indiv_sciddo_bed_files:
    # Optionally if you want merging this is the command
    # "bedtools sort | bedtools merge -d {params.merging_distance} | "
    input:
        expand(
            "results/sciddo/chromosomewise_chromHMM_model_{{sciddo_chromHMM_model}}/dump-bed/{{from_state}}-to-{{to_state}}/{{sciddo_group1}}-vs-{{sciddo_group2}}-{chromosome}.bed",
            chromosome = config["sciddo"]["chromatin_dynamics"]["chromosomes_to_examine"],
        ),
    output:
        sciddo_dynamics_merged_bed_out="results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/merged-bed/{from_state}-to-{to_state}-{sciddo_group1}-vs-{sciddo_group2}.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | "
        "bedtools sort | "
        "awk 'BEGIN{{OFS=\"\\t\"}} {{str = sprintf(\"%s:%s-%s\", $1, $2, $3)}} {{print $1, $2, $3, str}}' > {output.sciddo_dynamics_merged_bed_out}"

