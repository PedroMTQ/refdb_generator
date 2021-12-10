Provides automated creation of multiple databases:
- Rhea HMMs **implementation not finished**
- Reactome HMMs **implementation not finished**
- Swissprot diamond database 
- Trembl diamond database 
- Enzyme EC HMMs 
- BIGG reactions HMMs **implementation not finished**
- BIGG genes diamond database 


BIGG is generally used for metabolic modelling and these HMMs can thus be linked to metabolic models and metabolic model generation tools

All databases contain cross-linking metadata which is in a format compatible to Mantis - https://github.com/PedroMTQ/mantis

Arguments:
- `-db` or `--database`,to choose whether to cluster with `rhea`, `reactome`, `swissprot`, `trembl`, `ec`, `bigg_reactions` and `bigg_genes`
- `-o` or `--output_folder`,to specify directory to save database in
- `-ms` or `--min_seqs`, minimum sequences per HMM. Default is 10
- `-rf` or `--remove_files`, to remove files from previous runs



The `refdb_generator.yml` provides a conda environment recipe with all the required packages.

To run this tool:

1. `git clone git@github.com:PedroMTQ/hmm_updater.git`  
2. Go to cloned mantis folder and run `conda env create -f reference_generator.yml`
3. Run `conda activate reference_generator`
5. Run `python Reference_Generator.py` with the required arguments



Software used:
- HMMER to build and press HMMs
- Muscle for MSAs under 500 sequences
- Clustal Omega for MSAs above 500 sequences
- Diamond to build Diamond databases

Disclaimer:
The author does not own any of the third-party tools or data.
Please cite the respective tools and databases.