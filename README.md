Provides automated creation of multiple databases:
- Rhea HMMs **implementation not finished**
- Reactome diamond database
- Swissprot diamond database
- Trembl diamond database
- Enzyme EC HMMs
- BIGG genes diamond database


BIGG is generally used for metabolic modelling and these HMMs can thus be linked to metabolic models and metabolic model generation tools

All databases contain cross-linking metadata which is in a format compatible to Mantis - https://github.com/PedroMTQ/mantis

Arguments:
- `-db` or `--database`,to choose whether to cluster with `rhea`, `reactome`, `swissprot`, `trembl`, `ec`, and `bigg_genes`
- `-o` or `--output_folder`,to specify directory to save database in
- `-ms` or `--min_seqs`, minimum sequences per HMM. Default is 10
- `-rf` or `--remove_files`, to remove files from previous runs



The `refdb_generator.yml` provides a conda environment recipe with all the required packages.

To run this tool:

1. `git clone git@github.com:PedroMTQ/hmm_updater.git`  
2. Go to cloned mantis folder and run `conda env create -f reference_generator.yml`
3. Run `conda activate reference_generator`
5. Run `python Reference_Generator.py` with the required arguments

### How are these references created?

The reference HMMs were created by clustering protein sequences based on a certain database ID (e.g., for enzyme ECs, we created HMMs from the protein sequences that were annotated as having the enzyme EC function), whereas the diamond databases were created to be able to do a sequence homology search and then inferring all the associated functions of the match to the unknown sequence.

Specifically, the reference HMMs were created in the following manner:
\begin{enumerate}
    \item designate the clustering ID type (e.g. Rhea reactions)
    \item extract all protein sequences from the database - extracted directly or predicted with Prodigal \cite{hyatt_prodigal_2010}
    \item create fasta files where each file contains the protein sequences associated with a certain ID (e.g., all protein sequences associated with a certain Rhea reaction)
    \item Cluster the protein sequences in each fasta using MMSeqs2 \cite{mmseqs2}
    \item create new fasta files based on clustering results
    \item run a multiple sequence alignment on previous fasta files (MUSCLE \cite{muscle} for fasta files with 500 or less sequences, and Clustal Omega \cite{clustalomega} got files with more than 500 sequences)
    \item create HMMs using HMMER \cite{roberts_eddy_hmmer_nodate}
    \item index HMMs and create a corresponding metadata.tsv file
\end{enumerate}

Diamond databases were created by extracting the SwissProt/Trembl sequences, and then creating an associated metadata.tsv file.


### Software used
- HMMER to build and press HMMs
- Muscle for MSAs under 500 sequences
- Clustal Omega for MSAs above 500 sequences
- Diamond to build Diamond databases

Disclaimer:
The author does not own any of the third-party tools or data.
Please cite the respective tools and databases.
