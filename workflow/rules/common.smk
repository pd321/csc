import os
import datetime
import pandas as pd
from itertools import combinations, product, zip_longest
from snakemake.utils import validate, available_cpu_count
import time

report: "../report/workflow.rst"

configfile: 'config/config.yaml'
# validate(config, schema="schemas/config_schema.yaml")

chromhmm_metadata_file = "config/chromHMM_metadata.tsv"
metadata_df = pd.read_csv(chromhmm_metadata_file,sep="\t", index_col=0)

cell_types = config['chromHMM']['cell_types'] if "cell_types" in config['chromHMM'] else metadata_df.index.tolist()

# Setup vars
config_threads = int(config["threads"])
run_name = os.path.basename(os.getcwd())

sciddo_sample_types = config['sciddo']['sample_types']

sciddo_sample_combinations = list(
    itertools.combinations(iterable=sciddo_sample_types, r=2)
)

run_modes = config['run_mode'].split(',') if config['run_mode'] else []

def get_output():
	outfiles = expand("results/chromHMM/learn_model/{num_of_states}/model_{num_of_states}.txt", num_of_states = config['chromHMM']['num_of_states'])
	if "sciddo" in run_modes:
		sciddo_sample_and_model_combinations = list(
			itertools.product(
				config["sciddo"]["do_sciddo_for_these_models"], sciddo_sample_combinations
			)
		)
		outfiles += [
			f"results/sciddo/chromHMM_model_{model}/dump/{sample1}-vs-{sample2}.bed"
			for model, (sample1, sample2) in sciddo_sample_and_model_combinations
		]
	
	if "sagaconf" in run_modes:
		outfiles += expand("results/sagaCONF/run_results/{num_of_states}/{base_cell_type}-vs-{verif_cell_type}/confident_segments_dense.bed", num_of_states = config['sagaconf']['do_sagaconf_for_these_models'], base_cell_type = "CelltypeA_rep1", verif_cell_type = "CelltypeA_rep2")

	if "chromatin_dynamics" in run_modes:
		sciddo_chromatin_dynamics_combinations = list(
			itertools.product(
				config["sciddo"]["do_sciddo_for_these_models"],
				[config["sciddo"]["chromatin_dynamics"]["from_state"]],
				[config["sciddo"]["chromatin_dynamics"]["to_state"]],
				sciddo_sample_combinations,
			)
		)

		outfiles += [
			f"results/sciddo/chromosomewise_chromHMM_model_{sciddo_chromHMM_model}/merged-bed/{from_state}-to-{to_state}-{sciddo_group1}-vs-{sciddo_group2}.bed"
			for sciddo_chromHMM_model, from_state, to_state, (
				sciddo_group1,
				sciddo_group2,
			) in sciddo_chromatin_dynamics_combinations
		]
	return outfiles