# This will take all the separate files specified as input below and put them into the sciddo h5 file
# This file can then be used as needed for further analysis
rule convert_states_and_score:
    input:
        model_emissions="results/chromHMM/learn_model/{sciddo_chromHMM_model}/model_{sciddo_chromHMM_model}.txt",
        chrom_sizes=config["genome_sizes"],
        state_labels="config/sciddo_state_labels_for_state_{sciddo_chromHMM_model}.tsv",
        state_colours="config/sciddo_state_colours_for_state_{sciddo_chromHMM_model}.tsv",
        design_matrix="config/sciddo_design_matrix_for_state_{sciddo_chromHMM_model}.tsv",
    output:
        sciddo_convert_out="results/sciddo/chromHMM_model_{sciddo_chromHMM_model}/convert/sciddo_convert.h5",
    conda:
        "../envs/sciddo.yaml"
    threads: config["threads"]
    shell:
        "sciddo.py "
        "--no-conf-dump "
        "--workers {threads} "
        "convert "
        "--state-seg results/chromHMM/post_pca_subset/learn_model/{wildcards.sciddo_chromHMM_model}/ "
        "--chrom-sizes {input.chrom_sizes} "
        "--state-labels {input.state_labels} "
        "--state-colors {input.state_colours} "
        "--design-matrix {input.design_matrix} "
        "--model-emissions {input.model_emissions} "
        "--output {output.sciddo_convert_out} "
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

rule scan:
    input:
        sciddo_convert_out=rules.convert_states_and_score.output.sciddo_convert_out,
    output:
        sciddo_scan_out="results/sciddo/chromHMM_model_{sciddo_chromHMM_model}/scan/{sciddo_group1}-vs-{sciddo_group2}.h5",
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

rule dump:
    input:
        sciddo_scan_out=rules.scan.output.sciddo_scan_out,
    output:
        sciddo_dump_out="results/sciddo/chromHMM_model_{sciddo_chromHMM_model}/dump/{sciddo_group1}-vs-{sciddo_group2}.tsv",
    conda:
        "../envs/sciddo.yaml"
    threads: config["threads"]
    shell:
        "sciddo.py "
        "--no-conf-dump "
        "dump "
        "--data-file {input.sciddo_scan_out} "
        "--output {output.sciddo_dump_out} "
        "--data-type segments"

rule dump_to_bed:
    input:
        sciddo_dump_out=rules.dump.output.sciddo_dump_out,
    output:
        sciddo_dump_bed_out="results/sciddo/chromHMM_model_{sciddo_chromHMM_model}/dump/{sciddo_group1}-vs-{sciddo_group2}.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "tail -n +2 {input.sciddo_dump_out} | "
        "cut -f1,2,3,4,5 | "
        "bedtools sort > {output.sciddo_dump_bed_out}"
