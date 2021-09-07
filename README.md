Provides automated update of UniProt based HMMs, using two types of IDs as clustering, i.e., Rhea and ECs
Please set the directory paths where the HMMs should be created


The Rhea HMMs also contain cross-linking metadata which is in a format compatible to Mantis - https://github.com/PedroMTQ/mantis

Arguments:
- `-db` or `--database`,to choose whether to cluster with `ec` or  `rhea`
- `-o` or `--output_folder`,to specify directory to save HMMs in
- `-ms` or `--min_seqs`, minimum sequences per HMM. Default is 10
- `-rf` or `--remove_files`, to remove files from previous runs


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