rule move_posteriors:
	input:
		posterior = expand("results/chromHMM/learn_model/{{num_of_states}}/POSTERIOR/{{cell_types}}_{{num_of_states}}_{chromosome}_posterior.txt", chromosome = config['chromHMM']['chromosomes'])
	output:
		shifted_posteriors = expand("results/sagaCONF/posteriors/{{num_of_states}}/{{cell_types}}/{{cell_types}}_{{num_of_states}}_{chromosome}_posterior.txt", chromosome = config['chromHMM']['chromosomes'])
	shell:
		'mv {input} results/sagaCONF/posteriors/{wildcards.num_of_states}/{wildcards.cell_types}/'

rule sagaconf_parser:
	input:
		posteriors = rules.move_posteriors.output.shifted_posteriors
	output:
		merged_posterior = "results/sagaCONF/parsed_posteriors/{num_of_states}/{cell_types}/parsed_posterior.bed"
	conda:
		"../envs/sagaconf.yaml"
	shell:
		'python assets/bin/sagaconf/SAGAconf_parser.py --saga chmm results/sagaCONF/posteriors/{wildcards.num_of_states}/{wildcards.cell_types}/ 200 results/sagaCONF/parsed_posteriors/{wildcards.num_of_states}/{wildcards.cell_types}/'

rule sagaconf:
	input:
		merged_posterior_base = "results/sagaCONF/parsed_posteriors/{num_of_states}/{base_cell_type}/parsed_posterior.bed",
		merged_posterior_verif = "results/sagaCONF/parsed_posteriors/{num_of_states}/{verif_cell_type}/parsed_posterior.bed"
	output:
		confident_segments = "results/sagaCONF/run_results/{num_of_states}/{base_cell_type}-vs-{verif_cell_type}/confident_segments_dense.bed"
	conda:
		"../envs/sagaconf.yaml"
	shell:
		'python assets/bin/sagaconf/SAGAconf.py {input.merged_posterior_base} {input.merged_posterior_verif} results/sagaCONF/run_results/{wildcards.num_of_states}/{wildcards.base_cell_type}-vs-{wildcards.verif_cell_type}/'