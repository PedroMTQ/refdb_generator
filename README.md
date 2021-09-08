Provides automated update of UniProt based HMMs, using several types of IDs as clustering, i.e., Rhea, ECs, Reactome and BIGG

Please set the directory paths where the HMMs should be created

BIGG is generally used for metabolic modelling and these HMMs can thus be linked to metabolic models and metabolic model generation tools

All references contain cross-linking metadata which is in a format compatible to Mantis - https://github.com/PedroMTQ/mantis

Arguments:
- `-db` or `--database`,to choose whether to cluster with `ec`, `rhea`, `reactome`,`bigg_reaction`, and `bigg_genes`
- `-o` or `--output_folder`,to specify directory to save output in
- `-ms` or `--min_seqs`, minimum sequences per HMM. Default is 10
- `-rf` or `--remove_files`, to remove files from previous runs

The `bigg_genes` database will be sequence specific and therefore we will generate a diamond database based on all the sequences. Each sequence will have an associated metadata (bigg genes reactions and these reactions metadata).

The other databases will be the result of MSA clustering and are then converted into an HMM, as such the corresponding metadata will correspond to the whole HMM.

Keep in mind the multiple sequence alignment (MSA) might take a few days to run for some specific alignments.

The `hmm_update.yml` provides a conda environment recipe with all the required packages.

To run this tool:

1. `git clone git@github.com:PedroMTQ/hmm_updater.git`  
2. Go to cloned mantis folder and run `conda env create -f hmm_update.yml`
3. Run `conda activate hmm_update`
5. Run `python HMM_Updater.py` with the required arguments



Software used:
- HMMER to build and press HMMs
- Muscle for MSAs under 500 sequences
- Clustal Omega for MSAs above 500 sequences

Disclaimer:
The authors do not own any of the third-party tools or data.
Please cite the respective tools and databases.