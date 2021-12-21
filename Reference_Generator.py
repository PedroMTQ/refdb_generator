import urllib.request as request
import shutil
from contextlib import closing
import gzip
import re
import os
import subprocess
from pathlib import Path
from multiprocessing import Process, current_process, cpu_count, Manager
from time import sleep
import os
import argparse
from sys import argv
from Web_Connector import Web_Connector
import requests
from gzip import open as gzip_open

__author__ = "Pedro Queirós"
__status__ = "Production"
__credits__ = ['Pedro Queirós']
SPLITTER='/'



class Reference_Generator():
    def __init__(self,work_dir,remove_files=False,min_seqs=10,number_cores=None):
        self.manager = Manager()
        self.queue = self.manager.list()
        self.remove_files=remove_files
        self.min_seqs=min_seqs
        if number_cores:
            self.worker_count=int(number_cores)
        else:
            self.worker_count=self.check_environment_cores()
        if not work_dir.startswith(SPLITTER): work_dir=f'{SPLITTER}{work_dir}'
        if not work_dir.endswith(SPLITTER): work_dir=f'{work_dir}{SPLITTER}'
        self.work_dir=work_dir
        self.fasta_dir = f'{self.work_dir}fastas{SPLITTER}'
        self.aln_dir = f'{self.work_dir}aln{SPLITTER}'
        self.hmm_dir = f'{self.work_dir}hmm{SPLITTER}'
        self.cluster_dir = f'{self.work_dir}cluster{SPLITTER}'
        self.tmp_dir = f'{self.work_dir}tmp{SPLITTER}'
        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)
        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir,self.cluster_dir,self.tmp_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)



    def run_command(self,command,shell=False):
        process = subprocess.run(command, shell=shell)
        return process

    def run_prodigal(self):
        fasta_folder=f'{self.fasta_dir}{SPLITTER}'
        fna_list=[i for i in os.listdir(fasta_folder) if i.endswith('.fna')]
        for fna in fna_list:
            faa=fna.replace('.fna','.faa_pro')
            prodigal_command=f'prodigal -i {fasta_folder}{fna} -a {fasta_folder}{faa}'
            self.run_command(prodigal_command,shell=True)

    def read_mmseqs_cluster_tsv(self,cluster_tsv_path):
        clusters={}
        with open(cluster_tsv_path) as file:
            for line in file:
                line=line.strip('\n')
                representative_seq,other_seq=line.split('\t')
                if representative_seq not in clusters: clusters[representative_seq]=set()
                clusters[representative_seq].add(other_seq)
        res={}
        c=0
        for representative_seq in clusters:
            res[c]={representative_seq}
            res[c].update(clusters[representative_seq])
            c+=1
        return res

    def generate_clusters_fasta(self,cluster,clusters,cluster_seqs):
        for c in clusters:
            cluster_fasta_path = f'{self.fasta_dir}{cluster}_{c}.faa'
            with open(cluster_fasta_path,'a+') as file:
                for seq_id in clusters[c]:
                    seq=cluster_seqs[seq_id]
                    seq_line=f'>{seq_id}\n{seq}\n'
                    file.write(seq_line)


    def split_fastas_by_clusters(self):
        fasta_folder=f'{self.fasta_dir}'.rstrip(SPLITTER)
        shutil.move(fasta_folder,f'{fasta_folder}_raw_fastas')
        Path(fasta_folder).mkdir(parents=True, exist_ok=True)
        for cluster in os.listdir(self.cluster_dir):
            cluster_tsv=f'{self.cluster_dir}{cluster}{SPLITTER}{cluster}.clu_cluster.tsv'
            cluster_fasta=f'{self.cluster_dir}{cluster}{SPLITTER}{cluster}.clu_all_seqs.fasta'
            clusters=self.read_mmseqs_cluster_tsv(cluster_tsv)
            cluster_seqs = self.read_protein_fasta(cluster_fasta)
            self.generate_clusters_fasta(cluster,clusters,cluster_seqs)


    def run_mmseqs(self):
        faa_list=[i for i in os.listdir(self.fasta_dir) if i.endswith('.faa')]
        for faa in faa_list:
            current_cluster_dir = f'{self.cluster_dir}{faa}{SPLITTER}'.replace('.faa','')
            Path(current_cluster_dir).mkdir(parents=True, exist_ok=True)

            db=faa.replace('.faa','.db')
            cluster=faa.replace('.faa','.clu')
            tsv=faa.replace('.faa','.tsv')
            #mmseqs_createdb=f'mmseqs createdb {self.fasta_dir}{faa} {current_cluster_dir}{db}'
            #self.run_command(mmseqs_createdb,shell=True)
            #mmseqs_cluster=f'mmseqs cluster {current_cluster_dir}{db} {current_cluster_dir}{cluster} {self.tmp_dir}'
            #self.run_command(mmseqs_cluster,shell=True)
            #mmseqs_tsv=f'mmseqs createtsv {current_cluster_dir}{db} {current_cluster_dir}{db} {current_cluster_dir}{cluster} {current_cluster_dir}{tsv} {self.tmp_dir}'
            #self.run_command(mmseqs_tsv,shell=True)

            mmseqs_tsv=f'mmseqs easy-cluster {self.fasta_dir}{faa}  {current_cluster_dir}{cluster} {self.tmp_dir}'
            self.run_command(mmseqs_tsv,shell=True)



    def download_diamond(self):
        diamond_url = 'http://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz'
        archive_path = f'{self.work_dir}diamond-linux64.tar.gz'
        with requests.get(diamond_url, stream=True) as r:
            with open(archive_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        shutil.unpack_archive(archive_path, extract_dir=self.work_dir)
        os.remove(archive_path)

    def processes_handler(self, target_worker_function, add_sentinels=True):
        '''
        this will first generate one process per worker, then we add sentinels to the end of the list which will basically tell us when the queue is empty
        if we need to add new work (e.g. when doing taxa annotation) we just add the new work to the start of the list
        '''
        # os.getpid to add the master_pid
        if len(self.queue)<self.worker_count: worker_count=len(self.queue)
        else: worker_count=self.worker_count
        processes = [Process(target=target_worker_function, args=(self.queue, os.getpid(),)) for _ in range(worker_count)]
        # adding sentinel record since queue can be signaled as empty when its really not
        if add_sentinels:
            for _ in range(worker_count):   self.queue.append(None)
        for process in processes:
            process.start()
        # we could manage the processes memory here with a while cycle
        for process in processes:
            process.join()
            # exitcode 0 for sucessful exists
            if process.exitcode != 0:
                sleep(5)
                print('Ran into an issue, check the log for details. Exitting!')
                os._exit(1)



    #######
    def get_slurm_value(self,wanted_val, regex_pattern):
        res = None
        slurm_job_id = os.environ.get('SLURM_JOBID')
        if slurm_job_id:
            process = subprocess.run('sacct -j ' + str(slurm_job_id) + ' -o ' + wanted_val, shell=True,stdout=subprocess.PIPE)
            wanted = re.search(regex_pattern, str(process.stdout))
            if wanted: res = wanted.group()
        return res

    def check_environment_cores(self):
        res = self.get_slurm_value('AllocCPUS', re.compile('\d+'))
        if res:
            if int(res):
                print('Cores allocated by slurm:', res)
                return int(res)
            else:
                res = cpu_count()
                print('Cores allocated:', res)
                return int(res)
        else:
            res = cpu_count()
            print('Cores allocated:', res)
            return int(res)

    #low memory footprint_version
    def read_protein_fasta_generator(self,protein_fasta_path):
        query=None
        seq=[]
        if not os.path.exists(protein_fasta_path): return res
        with open(protein_fasta_path, 'r') as file:
            line = file.readline()
            while line:
                if line.startswith('>'):
                    if query:
                        if '|' in query:
                            query = query.split('|')[1]
                        if protein_fasta_path.endswith('_pro'):
                            query = query.split()[0]
                        seq=''.join(seq).upper()
                        if seq:
                            yield query,seq
                        seq=[]
                    query=line.replace('>','').strip()
                else:
                    seq.append(line.strip())
                line = file.readline()
            if query:
                if '|' in query:
                    query=query.split('|')[1]
                    if protein_fasta_path.endswith('_pro'):
                        query = query.split()[0]
                seq = ''.join(seq).upper()
                if seq:
                    yield query, seq

    #low memory footprint_version
    def read_protein_fasta(self,protein_fasta_path):
        query=None
        seq=[]
        res={}
        if not os.path.exists(protein_fasta_path): return res
        with open(protein_fasta_path, 'r') as file:
            line = file.readline()
            while line:
                if line.startswith('>'):
                    if query:
                        if '|' in query:
                            query = query.split('|')[1]
                        if protein_fasta_path.endswith('_pro'):
                            query = query.split()[0]
                        res[query]= ''.join(seq).upper()
                        seq=[]
                    query=line.replace('>','').strip()
                else:
                    seq.append(line.strip())
                line = file.readline()
            if query:
                if '|' in query:
                    query=query.split('|')[1]
                if protein_fasta_path.endswith('_pro'):
                    query = query.split()[0]
                res[query] = ''.join(seq).upper()
        return res

    def get_seqs_count(self,target_sample):
        total_seqs = 0
        with open(target_sample) as file:
            line = file.readline()
            while line:
                if line[0] == '>': total_seqs += 1
                line = file.readline()
        return total_seqs

    def download_file_ftp(self,url, file_path):
        with closing(request.urlopen(url)) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r, f)

    def unpack_gz(self,gz_file, output_file):
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    def msa(self,fasta_file):
        outfile = fasta_file.replace('.faa', '.aln')
        outfile_path = f'{self.aln_dir}{os.path.basename(outfile)}'
        print(fasta_file,outfile_path)
        seqs_count = self.get_seqs_count(fasta_file)

        if seqs_count >= self.min_seqs:
            if seqs_count == 1:
                shutil.copyfile(fasta_file, outfile_path)
            else:
                # muscle for small-mid aligns
                command=None
                if seqs_count <= 500:
                    command = f'muscle -in {fasta_file} -clwout {outfile_path} -clwstrict'
                # clustal omega for larger
                else:
                    command = f'clustalo -i {fasta_file} -o {outfile_path} --seqtype Protein --outfmt clu'
                if command:
                    print('Running command:',command)
                    subprocess.run(command, shell=True)

    def msa_worker_function(self,queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            self.msa(record)

    def launch_fastas_msa(self):
        completed_msa = os.listdir(self.aln_dir)
        for file in os.listdir(self.fasta_dir):
            if file.endswith('.faa'):
                file_aln = file.replace('.faa', '.aln')
                if file_aln not in completed_msa:
                    self.queue.append(f'{self.fasta_dir}{file}')

        self.processes_handler(self.msa_worker_function)

    def hmm_builder(self,aln_file):
        software = 'hmmbuild'
        outfile = aln_file.replace('.aln', '.hmm')
        outfile_path = f'{self.hmm_dir}{os.path.basename(outfile)}'
        command = f'{software} {outfile_path} {aln_file}'
        print('Running command:', command)
        subprocess.run(command, shell=True)

    def hmm_builder_worker_function(self,queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            self.hmm_builder(record)

    def launch_aln_hmmer(self):
        completed_hmm = os.listdir(self.hmm_dir)
        for file in os.listdir(self.aln_dir):
            if file.endswith('.aln'):
                file_aln = file.replace('.aln', '.hmm')
                if file_aln not in completed_hmm:
                    self.queue.append(f'{self.aln_dir}{file}')
        self.processes_handler(self.hmm_builder_worker_function)

    def concat_files(self,output_file, list_file_paths):
        print('Concatenating files into ', output_file)
        with open(output_file, 'wb') as wfd:
            for f in list_file_paths:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                # forcing disk write
                wfd.flush()
                os.fsync(wfd.fileno())

    def merge_profiles(self, output_file):
        if not os.path.exists(output_file):
            print('Merging profiles in ', self.hmm_dir)
            profiles = [self.hmm_dir + SPLITTER + i for i in os.listdir(self.hmm_dir) if i.lower().endswith('.hmm')]
            self.concat_files(output_file, profiles)
            print('Pressing profile', output_file)
            self.hmm_presser(output_file)

    def hmm_presser(self,hmm_file):
        software = 'hmmpress'
        command = f'{software} {hmm_file}'
        print('Running command:', command)
        subprocess.run(command, shell=True)

    def parse_bigg(self,rhea2bigg_path,wanted_dbs=[]):
        print('Parsing BIGG metadata')
        with open(rhea2bigg_path) as file:
            file.readline()
            line = file.readline()
            res = {}
            not_added = set()
            while line:
                line = line.strip('\n')
                bigg_id, name, reaction_string, model_list, database_links, old_bigg_ids = line.split('\t')
                if database_links:
                    database_links = database_links.split(';')
                    database_links = [i.strip() for i in database_links]
                    for db_link in database_links:
                        db, db_id = db_link.split(': ')
                        if db == 'RHEA' and 'rhea' in wanted_dbs and wanted_dbs:
                            db_id = db_id.split('/')[-1]
                            if db_id not in res: res[db_id] = set()
                            res[db_id].add(bigg_id)
                        elif db == 'Reactome Reaction' and 'reactome' in wanted_dbs and wanted_dbs:
                            db_id = db_id.split('/')[-1]
                            if db_id not in res: res[db_id] = set()
                            res[db_id].add(bigg_id)
                        elif db == 'EC Number' and 'ec' in wanted_dbs and wanted_dbs:
                            db_id = db_id.split('/')[-1]
                            if not db_id.endswith('-'):
                                if db_id not in res: res[db_id] = set()
                                res[db_id].add(bigg_id)
                        elif not wanted_dbs:
                            if db == 'RHEA': db='rhea'
                            elif db == 'Reactome Reaction': db='reactome'
                            elif db == 'EC Number': db='enzyme_ec'
                            elif db == 'MetaNetX (MNX) Equation': db='metanetx'
                            elif db == 'SEED Reaction': db='seed'
                            elif db == 'BioCyc': db='biocyc_reaction'
                            elif db == 'KEGG Reaction': db='kegg_reaction'
                            else:
                                print(db)
                            db_id = db_id.split('/')[-1]
                            if bigg_id not in res: res[bigg_id]={}
                            if db not in res[bigg_id]: res[bigg_id][db]=set()
                            res[bigg_id][db].add(db_id)
                        else:
                            not_added.add(db)
                line = file.readline()
        return res


#MSAs need to be checked for quality
class Reference_Generator_Uniprot(Reference_Generator):
    def __init__(self,work_dir,remove_files,min_seqs,number_cores,db):
        Reference_Generator.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs,number_cores=number_cores)
        self.db=db
        self.workflow_function()

    def parse_seq(self,file):
        res=[]
        line=file.readline()
        line_type, line_value = line[0:2], line[5:].strip().strip('\n')
        res.append(line_value)
        while not line.startswith('//'):
            line_type, line_value = line[0:2], line[5:].strip().strip('\n')
            res.append(line_value)
            line = file.readline()
        res=''.join(res)
        res=res.replace(' ','')
        res=res.strip()
        return res

    def parse_function(self,file):
        res = []
        line = file.readline()
        line_type, line_value = line[0:2], line[5:].strip().strip('\n')
        start_recording=False
        while line_type == 'CC':
            if line_value.startswith('-!- FUNCTION:'):
                temp=[]
                temp.append(line_value)
                start_recording=True
            elif start_recording and not line_value.startswith('-!-'):
                temp.append(line_value)
            elif start_recording and line_value.startswith('-!-'):
                temp=' '.join(temp)
                temp=temp.replace('-!- FUNCTION:','')
                temp=temp.split('{ECO:')[0]
                temp=temp.strip()
                isoform_search=re.search('\[Isoform \d+\]:',temp)
                if isoform_search:
                    isoform_search=isoform_search.group()
                    temp=temp.replace(isoform_search,'')
                pubmed_search=re.findall('PubMed:\d+,?\s?',temp)
                if pubmed_search:
                    for i in pubmed_search:
                        temp=temp.replace(i,'')
                    temp=temp.replace('()','')
                    temp=temp.replace(' .','.')
                temp=temp.strip()
                res.append(temp)
                start_recording=False
            line_type, line_value = line[0:2], line[5:].strip().strip('\n')
            line = file.readline()
        return res

    def yield_entries(self,uniprot_data):
        seq_id, db_type, taxon_id, sequence = None, None, None,None
        res={}
        wanted_db_types=[
            'eggnog',
            'go',
            'embl',
            'kegg_gene',
            'ncbi_gene',
            'pfam',
            'interpro',
            'hamap',
            'refseq',
            'panther',
            'ensemblbacteria',
            'antibodypedia',
            'orthodb',
            'reactome',
            'genecards',
            'malacards',
            'prints',
            'tigrfams',
            'prosite',
            'string',
            'biogrid',
            'smart',
            'gene3d',
            'peptideatlas',
            'proteomicsdb',
            'phosphositeplus',
            'biocyc',
            'brenda',
            'pdb',
            'ensembl',
            'ensemblmetazoa',
            'ensemblplants',
            'expressionatlas',
            'genevisible',
            'ensemblfungi',
            'plantreactome',
            'swisslipids',
            'pathwaycommons',
            'chembl',
            'patric',
            'genewiki',
            'moondb',
                         ]
        with open(uniprot_data) as file:
            for line in file:
                db_type = None
                line=line.strip('\n')
                line_type, line_value = line[0:2], line[5:].strip()
                if line_type=='ID':
                    db_type=None
                    if seq_id:
                        if 'eggnog' in res:
                            for eggnog_id in res['eggnog']:
                                if re.search('arCOG\d+',eggnog_id):
                                    if 'arcog' not in res: res['arcog']=set()
                                    res['arcog'].add(eggnog_id)
                                elif re.search('COG\d+',eggnog_id):
                                    if 'cog' not in res: res['cog']=set()
                                    res['cog'].add(eggnog_id)
                        yield seq_id,taxon_id,res,sequence
                    res={}
                    seq_id,db_type,taxon_id,sequence = None,None,None,None
                    seq_id=line_value.split()[0]

                elif line_type == 'AC':
                    line_value=line_value.split(';')
                    db_type='uniprot'
                elif line_type == 'GN' and line_value.startswith('Name='):
                    db_type = 'uniprot_gene'
                    line_value = line_value.replace('Name=', '').split(';')[0]
                    line_value = line_value.split('{ECO:')[0]
                elif line_type == 'OX' and line_value.startswith('NCBI_TaxID='):
                    db_type = None
                    line_value=line_value.replace('NCBI_TaxID=','').strip()
                    line_value = line_value.split('{ECO:')[0]
                    line_value=line_value.strip(';').strip()
                    taxon_id=line_value

                elif line_type == 'DE':
                    db_type=None
                    if line_value.startswith('RecName') or line_value.startswith('AltName'):
                        db_type = 'description'
                        line_value=line_value.replace('RecName:','').strip()
                        line_value=line_value.replace('AltName:','').strip()
                        line_value=line_value.replace('Full=','').strip()
                        line_value=line_value.split('{ECO:')[0]
                        line_value=line_value.strip(';')
                    elif line_value.startswith('EC='):
                        db_type = 'enzyme_ec'
                        line_value=line_value.replace('EC=','').strip().split()[0]
                        line_value=line_value.strip(';')
                elif line_type=='CC':
                    db_type = 'description'
                    line_value=self.parse_function(file)
                elif line_type=='SQ':
                    db_type=None
                    sequence=self.parse_seq(file)
                elif line_type=='DR':
                    #database reference
                    db_type,line_value= line_value.split(';')[0:2]
                    db_type,line_value=db_type.strip(),line_value.strip()
                    db_type=db_type.lower()
                    if db_type=='go':
                        line_value=line_value.replace('GO:','')
                    elif db_type=='kegg': db_type='kegg_gene'
                    elif db_type=='geneid': db_type='ncbi_gene'
                    elif db_type=='tigrfams': db_type='tigrfam'
                    if db_type not in wanted_db_types: db_type=None
                if db_type:
                    if not isinstance(line_value,list): line_value=[line_value]
                    if db_type not in res: res[db_type]=set()
                    line_value=[i.strip() for i in line_value]
                    line_value=[i for i in line_value if i]
                    res[db_type].update(line_value)
        yield seq_id, taxon_id, res, sequence

    def output_uniprot_data_hmm(self,data_generator,db_type):
        for seq_id,taxon_id,seq_metadata,sequence in data_generator:
            if db_type in seq_metadata:
                for db_id in seq_metadata[db_type]:
                    current_fasta = f'{self.fasta_dir}{db_id}.faa'
                    with open(current_fasta, 'a+') as file:
                        line = f'>{seq_id}\n{sequence}'
                        file.write(f'{line}\n')
        self.run_mmseqs()
        self.split_fastas_by_clusters()

    def create_hmms(self,db_type):
        data_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
        compressed_data = f'{self.work_dir}uniprot_sprot.dat.gz'
        uncompressed_data = f'{self.work_dir}uniprot_sprot.dat'
        if not os.path.exists(compressed_data) or not os.path.exists(uncompressed_data):
            print(f'Downloading data:\n{data_url}')
            self.download_file_ftp(data_url, compressed_data)

        if not os.path.exists(uncompressed_data):
            print(f'Uncompressing data:\n{compressed_data}')
            self.unpack_gz(compressed_data, uncompressed_data)
        data_generator=self.yield_entries(uncompressed_data)
        self.output_uniprot_data_hmm(data_generator,db_type)
        #self.launch_fastas_msa()
        #self.launch_aln_hmmer()
        #hmm_file=f'{self.work_dir}uniprot_ec.hmm'
        #self.merge_profiles(output_file=hmm_file)
        #print(f'Finished generating {hmm_file}')

    def output_uniprot_data_taxa(self,data_generator):
        general_taxon=self.get_ncbi_domains()
        for seq_id,taxon_id,seq_metadata,sequence in data_generator:
            c+=1
            if not taxon_id:
                current_folder=f'{self.fasta_dir}uniprotG{SPLITTER}'
                Path(current_folder).mkdir(parents=True, exist_ok=True)
                current_fasta=f'{current_folder}uniprotG_merged.faa'
                current_metadata=f'{current_folder}metadata.tsv'
            else:
                current_folder=f'{self.fasta_dir}{taxon_id}{SPLITTER}'
                Path(current_folder).mkdir(parents=True, exist_ok=True)
                current_fasta=f'{current_folder}{taxon_id}_merged.faa'
                current_metadata=f'{current_folder}metadata.faa'
            with open(current_fasta, 'a+') as file:
                line = f'>{seq_id}\n{sequence}'
                file.write(f'{line}\n')

            with open(current_metadata, 'a+') as file:
                temp = [f'{seq_id}','|']
                for db in seq_metadata:
                    for db_id in seq_metadata[db]:
                        temp.append(f'{db}:{db_id}')
                line='\t'.join(temp)
                file.write(f'{line}\n')

    def output_uniprot_data_dmnd(self,data_generator,fasta_file,metadata_file):
        if os.path.exists(fasta_file):
            os.remove(fasta_file)
        if os.path.exists(metadata_file):
            os.remove(metadata_file)

        for seq_id,taxon_id,seq_metadata,sequence in data_generator:
            with open(fasta_file, 'a+') as file:
                line = f'>{seq_id}\n{sequence}'
                file.write(f'{line}\n')
            with open(metadata_file, 'a+') as file:
                temp = [f'{seq_id}', '|']
                for db in seq_metadata:
                    for db_id in seq_metadata[db]:
                        temp.append(f'{db}:{db_id}')
                line = '\t'.join(temp)
                file.write(f'{line}\n')

    def create_diamond_dbs(self):
        if self.db=='swissprot':
            data_url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
            compressed_data = f'{self.work_dir}uniprot_sprot.dat.gz'
            uncompressed_data = f'{self.work_dir}uniprot_sprot.dat'
        elif self.db=='trembl':
            data_url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz'
            compressed_data = f'{self.work_dir}uniprot_trembl.dat.gz'
            uncompressed_data = f'{self.work_dir}uniprot_trembl.dat'

        if not os.path.exists(compressed_data) or not os.path.exists(uncompressed_data):
            print(f'Downloading data:\n{data_url}')
            self.download_file_ftp(data_url, compressed_data)

        if not os.path.exists(uncompressed_data):
            print(f'Uncompressing data:\n{compressed_data}')
            self.unpack_gz(compressed_data, uncompressed_data)

        data_generator=self.yield_entries(uncompressed_data)

        fasta_file = f'{self.fasta_dir}uniprot_merged.faa'
        metadata_file = f'{self.fasta_dir}metadata.tsv'
        self.output_uniprot_data_dmnd(data_generator,fasta_file,metadata_file)

        self.download_diamond()
        diamond_path = f'{self.work_dir}diamond'
        dmnd_file=fasta_file.replace('.faa','')
        dmnd_command = f'{diamond_path} makedb --in {fasta_file} -d {dmnd_file}'
        subprocess.run(dmnd_command.split())
        os.remove(diamond_path)
        print(f'Finished generating {dmnd_file}')

    def workflow_function(self):
        if self.db=='swissprot' or self.db=='trembl':
            self.create_diamond_dbs()
        elif self.db=='ec':
            self.create_hmms(db_type='enzyme_ec')

#MSAs need to be checked for quality
class Reference_Generator_Rhea(Reference_Generator):
    def __init__(self,work_dir,remove_files,min_seqs,number_cores):
        Reference_Generator.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs,number_cores=number_cores)
        self.workflow_function()
        self.write_metadata()


    def parse_rhea2uniprot(self,rhea2uniprot_path):
        print('Parsing rhea2uniprot')
        res = {}
        with open(rhea2uniprot_path) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n')
                if line:
                    rhea_id, direction, master_id, uniprot_id = line.split('\t')
                    if master_id not in res: res[master_id] = set()
                    res[master_id].add(uniprot_id)
                line = file.readline()
        return res



    def write_metadata(self):


        rhea2xrefs_file = f'{self.work_dir}rhea2xrefs.tsv'
        rhea2xrefs_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2xrefs.tsv'

        rhea2bigg_file=f'{self.work_dir}bigg_models_reactions.txt'
        rhea2bigg_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'


        if not os.path.exists(rhea2xrefs_file):
            self.download_file_ftp(rhea2xrefs_url, rhea2xrefs_file)

        if not os.path.exists(rhea2bigg_file):
            self.download_file_ftp(rhea2bigg_url, rhea2bigg_file)

        wanted_db='rhea'
        metadata_file = f'{self.work_dir}metadata.tsv'
        rhea2bigg=self.parse_bigg(rhea2bigg_file,wanted_dbs=[wanted_db])


        if not os.path.exists(metadata_file):
            print('Parsing rhea2xrefs')
            rhea2ids={}
            with open(rhea2xrefs_file) as file:
                line=file.readline()
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    if line:
                        rhea_id,direction,master_id,db_id,db_type=line.split('\t')
                        if master_id not in rhea2ids: rhea2ids[master_id]={}
                        if db_type=='EC': db_type='enzyme_ec'
                        elif db_type=='METACYC': db_type='biocyc_reaction'
                        elif db_type=='ECOCYC': db_type='biocyc_reaction'
                        elif db_type=='KEGG_REACTION': db_type='kegg_reaction'
                        elif db_type=='GO':
                            db_type='go'
                            db_id=db_id.strip('GO:')

                        elif db_type=='REACTOME': db_type=None
                        elif db_type=='MACIE': db_type=None
                        if db_type:
                            if db_type not in rhea2ids[master_id]: rhea2ids[master_id][db_type]=set()
                            rhea2ids[master_id][db_type].add(db_id)
                    line=file.readline()
            with open(metadata_file,'w+') as file:
                for main_id in rhea2ids:

                    line = [main_id,'|']
                    line.append(f'{wanted_db}:{main_id}')
                    if main_id in rhea2bigg:
                        for bigg_id in rhea2bigg[main_id]:
                            line.append(f'bigg:{bigg_id}')
                    for db_type in rhea2ids[main_id]:
                        for db_id in rhea2ids[main_id][db_type]:
                            line.append(f'{db_type}:{db_id}')
                    file.write('\t'.join(line)+'\n')
        os.remove(rhea2xrefs_file)
        os.remove(rhea2bigg_file)

    def fasta_writer(self,uniprot_mapping, uncompressed_uniprot_fastas):
        uniprot_seqs=self.read_protein_fasta_generator(uncompressed_uniprot_fastas)
        for uniprot_seq in uniprot_seqs:
            uniprot_id,sequence=uniprot_seq
            for main_id in uniprot_mapping:
                seq_ids=uniprot_mapping[main_id]
                if uniprot_id in seq_ids:
                    fasta_file = f'{self.fasta_dir}{main_id}.faa'
                    with open(fasta_file, 'a+') as file:
                        outline = f'>{uniprot_id}\n{sequence}\n'
                        file.write(outline)


    def workflow_function(self):

        uniprot_fastas_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
        compressed_uniprot_fastas = f'{self.work_dir}uniprot_file.xml.gz'
        uncompressed_uniprot_fastas = f'{self.work_dir}uniprot_file.xml'

        rhea2uiniprot_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv'
        rhea2uiniprot_file = f'{self.work_dir}rhea2uniprot.tsv'


        hmm_file=f'{self.work_dir}uniprot_rhea.hmm'



        if not os.path.exists(compressed_uniprot_fastas) or not os.path.exists(uncompressed_uniprot_fastas):
            self.download_file_ftp(uniprot_fastas_url, compressed_uniprot_fastas)


        if not os.path.exists(rhea2uiniprot_file):
            self.download_file_ftp(rhea2uiniprot_url, rhea2uiniprot_file)




        if not os.path.exists(uncompressed_uniprot_fastas):
            self.unpack_gz(compressed_uniprot_fastas, uncompressed_uniprot_fastas)

        if not os.listdir(self.fasta_dir):
            rhea_uniprot = self.parse_rhea2uniprot(rhea2uiniprot_file)
            self.fasta_writer(rhea_uniprot, uncompressed_uniprot_fastas)

        self.launch_fastas_msa()
        self.launch_aln_hmmer()
        self.merge_profiles(output_file=hmm_file)
        print(f'Finished generating {hmm_file}')

class Reference_Generator_Reactome(Reference_Generator):
    def __init__(self,work_dir,remove_files,min_seqs,number_cores):
        Reference_Generator.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs,number_cores=number_cores)
        self.workflow_function()


    def parse_reactome2uniprot(self,reactome2uniprot_path):
        print('Parsing reactome2uniprot')
        res = {}
        with open(reactome2uniprot_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if line:
                    line = line.split('\t')
                    uniprot_id, master_id = line[0], line[1]
                    if uniprot_id not in res: res[uniprot_id] = set()
                    res[uniprot_id].add(master_id)
                line = file.readline()
        return res

    def write_metadata(self,reactome_uniprot):
        bigg2refs_file=f'{self.work_dir}bigg_models_reactions.txt'
        bigg2refs_file_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'

        if not os.path.exists(bigg2refs_file):
            self.download_file_ftp(bigg2refs_file_url, bigg2refs_file)
        wanted_db='reactome'
        metadata_file = f'{self.work_dir}metadata.tsv'
        bigg_metadata=self.parse_bigg(bigg2refs_file,wanted_dbs=[wanted_db])
        if not os.path.exists(metadata_file):
            with open(metadata_file,'w+') as file:
                for uniprot_id in reactome_uniprot:
                    line = [uniprot_id,'|']
                    for reactome_id in reactome_uniprot[uniprot_id]:
                        line.append(f'reactome:{reactome_id}')
                        if reactome_id in bigg_metadata:
                            for db_id in bigg_metadata[reactome_id]:
                                line.append(f'bigg:{db_id}')
                    file.write('\t'.join(line)+'\n')
        os.remove(bigg2refs_file)

    def merge_faa(self,main_fasta):
        fasta_folder=f'{self.fasta_dir}{SPLITTER}'
        with open(main_fasta, 'a+') as file:
            for fasta in os.listdir(fasta_folder):
                all_sequences=self.read_protein_fasta_generator(fasta)
                for seq_id,protein_sequence in all_sequences:
                    outline = f'>{seq_id}\n{protein_sequence}\n'
                    file.write(outline)

    def fasta_writer(self,fasta_path,uniprot_mapping, uncompressed_uniprot_fastas):
        uniprot_seqs=self.read_protein_fasta_generator(uncompressed_uniprot_fastas)
        with open(fasta_path, 'w+') as file:
            for uniprot_id,sequence in uniprot_seqs:
                if uniprot_id in uniprot_mapping:
                    outline = f'>{uniprot_id}\n{sequence}\n'
                    file.write(outline)


    def workflow_function(self):

        uniprot_fastas_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
        compressed_uniprot_fastas = f'{self.work_dir}uniprot_file.faa.gz'
        uncompressed_uniprot_fastas = f'{self.work_dir}uniprot_file.faa'

        reactome2uniprot_url='https://reactome.org/download/current/UniProt2ReactomeReactions.txt'
        reactome2uniprot_file = f'{self.work_dir}UniProt2ReactomeReactions.txt'

        dmnd_file=f'{self.work_dir}reactome'
        fasta_path = f'{self.work_dir}reactome.faa'


        if not os.path.exists(compressed_uniprot_fastas) or not os.path.exists(uncompressed_uniprot_fastas):
            self.download_file_ftp(uniprot_fastas_url, compressed_uniprot_fastas)


        if not os.path.exists(reactome2uniprot_file):
            self.download_file_ftp(reactome2uniprot_url, reactome2uniprot_file)

        if not os.path.exists(uncompressed_uniprot_fastas):
            self.unpack_gz(compressed_uniprot_fastas, uncompressed_uniprot_fastas)

        if not os.listdir(self.fasta_dir):
            reactome_uniprot = self.parse_reactome2uniprot(reactome2uniprot_file)
            self.fasta_writer(fasta_path,reactome_uniprot, uncompressed_uniprot_fastas)

        diamond_path=f'{self.work_dir}diamond'

        if not os.path.exists(diamond_path):
            self.download_diamond()
        dmnd_command=f'{diamond_path} makedb --in {fasta_path} -d {dmnd_file}'
        subprocess.run(dmnd_command.split())
        self.write_metadata(reactome_uniprot)
        os.remove(diamond_path)
        print(f'Finished generating {dmnd_file}.dmnd')
        os.remove(reactome2uniprot_file)
        os.remove(compressed_uniprot_fastas)
        os.remove(uncompressed_uniprot_fastas)

class Reference_Generator_BIGG(Reference_Generator,Web_Connector):
    def __init__(self,work_dir,remove_files,min_seqs,number_cores):
        Reference_Generator.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs,number_cores=number_cores)
        Web_Connector.__init__(self)
        self.mp_results = self.manager.list()
        self.workflow_function()


    def get_all_models(self):
        res=set()
        models_url='http://bigg.ucsd.edu/api/v2/models'
        json_page=self.get_url_json(models_url)
        results=json_page['results']
        for i in results:
            bigg_id=i['bigg_id']
            res.add(bigg_id)
        return res

    def get_genes_model(self,model_id):
        res=set()
        models_url=f'http://bigg.ucsd.edu/api/v2/models/{model_id}/genes'
        json_page=self.get_url_json(models_url)
        results=json_page['results']
        for i in results:
            bigg_id=i['bigg_id']
            res.add(bigg_id)
        return res

    def get_gene_info(self,model_id,gene_id):
        models_url=f'http://bigg.ucsd.edu/api/v2/models/{model_id}/genes/{gene_id}'
        #print(f'Getting info for model {model_id} and gene {gene_id}')
        json_page=self.get_url_json(models_url)
        protein_sequence=json_page['protein_sequence']
        dna_sequence=json_page['dna_sequence']
        if protein_sequence or dna_sequence:
            reactions=json_page['reactions']
            reactions_bigg=set()
            for i in reactions:
                bigg_id=i['bigg_id']
                reactions_bigg.add(bigg_id)
            return [gene_id,protein_sequence,dna_sequence,reactions_bigg]

    def gene_info_worker_function(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            arg1, arg2 = record
            gene_info=self.get_gene_info(arg1, arg2)
            if gene_info:
                self.mp_results.append(gene_info)

    def launch_reaction_info_retrieval(self, model_id, genes_list):
        for gene_id in genes_list:
            self.queue.append([model_id, gene_id])
        self.processes_handler(self.gene_info_worker_function)
        while self.mp_results:
            yield self.mp_results.pop(0)

    def export_to_fasta(self,model_id,gene_id,protein_sequence,dna_sequence):
        fasta_path_aa = f'{self.fasta_dir}{model_id}.faa_pre'
        fasta_path_nt = f'{self.fasta_dir}{model_id}.fna'
        if protein_sequence:
            with open(fasta_path_aa, 'a+') as file:
                outline = f'>{gene_id}\n{protein_sequence}\n'
                file.write(outline)
        if dna_sequence:
            with open(fasta_path_nt, 'a+') as file:
                outline = f'>{gene_id}\n{dna_sequence}\n'
                file.write(outline)

    def write_metadata(self,bigg_metadata,metadata_file,sequences,prodigal_proteins):
        with open(metadata_file,'a+') as file:
            for seq_id in sequences:
                if seq_id in prodigal_proteins:
                    baseline_seq_id=seq_id.split('_')[0:-1]
                    baseline_seq_id='_'.join(baseline_seq_id)
                else: baseline_seq_id =seq_id
                line = [baseline_seq_id,'|']
                if baseline_seq_id in self.genes_reactions:
                    for reaction_id in self.genes_reactions[baseline_seq_id]:
                        line.append(f'bigg_reaction:{reaction_id}')
                        if reaction_id in bigg_metadata:
                            for db in bigg_metadata[reaction_id]:
                                for db_id in bigg_metadata[reaction_id][db]:
                                    line.append(f'{db}:{db_id}')
                    file.write('\t'.join(line)+'\n')
                else:
                    print(f'Seq missing {baseline_seq_id}')


    def merge_faa(self):
        '''
        this will merge the two fastas, one generated with the protein sequences extracted from bigg,
        the other with the protein sequences predicted with prodigal
        the resulting fasta will add protein sequences from bigg per gene, if not available it adds protein sequences from prodigal
        '''
        bigg_file=f'{self.work_dir}bigg_models_reactions.txt'
        bigg_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
        if not os.path.exists(bigg_file):
            self.download_file_ftp(bigg_url, bigg_file)
        metadata_file = f'{self.work_dir}metadata.tsv'
        bigg_metadata=self.parse_bigg(bigg_file)

        fasta_folder=f'{self.fasta_dir}{SPLITTER}'
        model_list=[i.replace('.fna','') for i in os.listdir(fasta_folder) if i.endswith('.fna')]
        for model_id in model_list:
            faa_pre=f'{fasta_folder}{model_id}.faa_pre'
            faa_prodigal=f'{fasta_folder}{model_id}.faa_pro'
            fna=f'{fasta_folder}{model_id}.fna'
            all_pre_proteins=self.read_protein_fasta(faa_pre)
            all_prodigal_proteins=self.read_protein_fasta(faa_prodigal)
            all_sequences=set(list(all_pre_proteins.keys())+list(all_prodigal_proteins.keys()))
            self.write_metadata(bigg_metadata,metadata_file,all_sequences,all_prodigal_proteins)
            for seq_id in all_sequences:
                fasta_path_aa = f'{self.fasta_dir}{SPLITTER}{model_id}.faa'
                protein_sequence=None
                #if bigg provides protein sequences we use them
                if seq_id in all_pre_proteins:
                    protein_sequence=all_pre_proteins[seq_id]
                #if bigg doesnt, we predict with prodigal
                elif seq_id not in all_pre_proteins and seq_id in all_prodigal_proteins:
                    protein_sequence=all_prodigal_proteins[seq_id]
                #if there is not protein and dna sequence, we just report it
                else:
                    print(f'Did not manage to export sequence {seq_id} for model {model_id}')
                if protein_sequence:
                    with open(fasta_path_aa, 'a+') as file:
                        outline = f'>{seq_id}\n{protein_sequence}\n'
                        file.write(outline)

    def export_bigg_faa(self):
        '''
        merges all the sequences into one big file
        '''
        fasta_folder=f'{self.fasta_dir}{SPLITTER}'
        faa_list=[i for i in os.listdir(fasta_folder) if i.endswith('.faa')]
        bigg_path = f'{self.work_dir}bigg.faa'
        with open(bigg_path, 'a+') as file:
            for faa in faa_list:
                faa_path=f'{fasta_folder}{faa}'
                all_seqs_generator=self.read_protein_fasta_generator(faa_path)
                for seq_id,protein_sequence in all_seqs_generator:
                        outline = f'>{seq_id}\n{protein_sequence}\n'
                        file.write(outline)

    def fasta_writer(self):
        self.genes_reactions={}
        for model_id in self.get_all_models():
            print(f'Getting info for model {model_id}')
            genes_list=self.get_genes_model(model_id)
            reactions_generator=self.launch_reaction_info_retrieval(model_id,genes_list)
            for gene_id, protein_sequence,dna_sequence, reactions_bigg in reactions_generator:
                if gene_id not in self.genes_reactions: self.genes_reactions[gene_id]=set()
                self.export_to_fasta(model_id,gene_id, protein_sequence,dna_sequence)
                for reaction_id in reactions_bigg:
                    self.genes_reactions[gene_id].add(reaction_id)
        self.run_prodigal()
        self.merge_faa()
        self.export_bigg_faa()


    def workflow_function(self):

        dmnd_file=f'{self.work_dir}bigg'
        fasta_path = f'{self.work_dir}bigg.faa'
        bigg_file=f'{self.work_dir}bigg_models_reactions.txt'


        self.fasta_writer()
        self.download_diamond()
        diamond_path=f'{self.work_dir}diamond'
        dmnd_command=f'{diamond_path} makedb --in {fasta_path} -d {dmnd_file}'
        subprocess.run(dmnd_command.split())
        print(f'Finished generating {dmnd_file}')

        os.remove(diamond_path)
        os.remove(bigg_file)


if __name__ == '__main__':
    print('Executing command:\n', ' '.join(argv))
    parser = argparse.ArgumentParser(description='This is a functional annotation reference generator tool\n'
                                     , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-db','--database', help='[required]\tClustering ID',choices=['ec', 'rhea','reactome','bigg_genes','swissprot','trembl'])
    parser.add_argument('-o', '--output_folder', help='[required]\tDirectory to save HMMs in')
    parser.add_argument('-c', '--number_cores', help='[optional]\tNumber of cores to use')
    parser.add_argument('-ms', '--min_seqs',help='[optional]\tMinimum sequences per HMM. Default is 10')
    parser.add_argument('-rf', '--remove_files', action='store_true',help='[optional]\tuse this to remove files from previous runs.')


    args = parser.parse_args()
    database = args.database
    output_folder = args.output_folder
    min_seqs = args.min_seqs
    number_cores = args.number_cores
    remove_files = args.remove_files
    if min_seqs:    min_seqs=int(min_seqs)
    else:           min_seqs=10
    if not output_folder:
        print('Missing output folder!')
    elif database=='rhea':
        updater = Reference_Generator_Rhea(output_folder, remove_files,min_seqs,number_cores)
    elif database=='reactome':
        updater=Reference_Generator_Reactome(output_folder,remove_files,min_seqs,number_cores)
    elif database=='swissprot':
        updater=Reference_Generator_Uniprot(output_folder,remove_files,min_seqs,number_cores,db='swissprot')
    elif database=='trembl':
        updater=Reference_Generator_Uniprot(output_folder,remove_files,min_seqs,number_cores,db='trembl')
    elif database=='ec':
        updater=Reference_Generator_Uniprot(output_folder,remove_files,min_seqs,number_cores,db='ec')
    elif database=='bigg_genes':
        updater=Reference_Generator_BIGG(output_folder,remove_files,min_seqs,number_cores)
    else:
        print('Command is not valid')
