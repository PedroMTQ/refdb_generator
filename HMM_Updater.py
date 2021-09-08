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


__author__ = "Pedro Queirós"
__status__ = "Production"
__credits__ = ['Pedro Queirós','Polina Novikova']
SPLITTER='/'



class HMM_Updater():
    def __init__(self,work_dir,remove_files=False,min_seqs=10):
        self.manager = Manager()
        self.queue = self.manager.list()
        self.remove_files=remove_files
        self.min_seqs=min_seqs
        worker_count = self.check_environment_cores()
        self.worker_count=worker_count
        self.work_dir=work_dir
        self.fasta_dir = f'{self.work_dir}{SPLITTER}fastas{SPLITTER}'
        self.aln_dir = f'{self.work_dir}{SPLITTER}aln{SPLITTER}'
        self.hmm_dir = f'{self.work_dir}{SPLITTER}hmm{SPLITTER}'


    def processes_handler(self, target_worker_function, add_sentinels=True):
        '''
        this will first generate one process per worker, then we add sentinels to the end of the list which will basically tell us when the queue is empty
        if we need to add new work (e.g. when doing taxa annotation) we just add the new work to the start of the list
        '''
        # os.getpid to add the master_pid
        processes = [Process(target=target_worker_function, args=(self.queue, os.getpid(),)) for _ in range(self.worker_count)]
        # adding sentinel record since queue can be signaled as empty when its really not
        if add_sentinels:
            for _ in range(self.worker_count):   self.queue.append(None)
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
        with open(protein_fasta_path, 'r') as file:
            line = file.readline()
            while line:
                if line.startswith('>'):
                    if query:
                        query = query.split('|')[1]
                        yield query,''.join(seq).upper()
                        seq=[]
                    query=line.replace('>','').strip()
                else:
                    seq.append(line.strip())
                line = file.readline()
            if query:
                query=query.split('|')[1]
                yield query, ''.join(seq).upper()


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
        seqs_count = self.get_seqs_count(fasta_file)
        if seqs_count >= self.min_seqs:
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


class HMM_Updater_Uniprot_EC(HMM_Updater):
    def __init__(self,work_dir,remove_files,min_seqs):
        HMM_Updater.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs)
        self.workflow_function()
        self.write_metadata()


    def uniprot_xml_parser(self,input_file_xml):
        dict_accession_seqs = dict()
        dict_ec_accession = dict()

        with open(input_file_xml, 'r') as file:
            line = file.readline()
            while line:
                line = line.strip('\n').strip()
                if line.startswith('<accession>'):
                    regex_pattern = re.compile('>([A-Z]|\d)+<')
                    accession = re.search(regex_pattern, line).group().strip('<>')
                elif line.startswith('<sequence'):
                    sequence = line.split('>')[1]
                    sequence = sequence.split('<')[0]
                    dict_accession_seqs[accession] = sequence

                elif line.startswith('<ecNumber'):
                    ec_number = line.split('>')[1]
                    ec_number = ec_number.split('<')[0]
                    if '-' not in ec_number and not re.search('[A-Za-z]',ec_number):
                        if ec_number not in dict_ec_accession:
                            dict_ec_accession[ec_number] = set()
                        dict_ec_accession[ec_number].add(accession)

                elif line.startswith('<dbReference') and 'type="EC"' in line:
                    db_ref = line.split('id="')[1]
                    db_ref = db_ref.split('"')[0]
                    if '-' not in db_ref and not re.search('[A-Za-z]',db_ref):
                        if db_ref not in dict_ec_accession:
                            dict_ec_accession[db_ref] = set()
                        dict_ec_accession[db_ref].add(accession)

                line = file.readline()

            return dict_accession_seqs, dict_ec_accession


    def fasta_writer(self,dict_accession_seqs, dict_ec_accession):
        for ec_number in dict_ec_accession:
            fasta_file = f'{self.fasta_dir}{ec_number}.faa'
            if len(dict_ec_accession[ec_number])>=10:
                with open(fasta_file, 'w+') as file:
                    for accession in dict_ec_accession[ec_number]:
                        if accession in dict_accession_seqs:
                            sequence = dict_accession_seqs[accession]
                            outline = f'>{accession}\n{sequence}\n'
                            file.write(outline)

    def write_metadata(self):
        bigg2refs_file=f'{self.work_dir}{SPLITTER}bigg_models_reactions.txt'
        bigg2refs_file_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'

        if not os.path.exists(bigg2refs_file):
            self.download_file_ftp(bigg2refs_file_url, bigg2refs_file)

        wanted_ec='ec'
        metadata_file = f'{self.work_dir}{SPLITTER}uniprot_{wanted_ec}.tsv'
        bigg_metadata=self.parse_bigg(bigg2refs_file,wanted_dbs=['ec'])

        if not os.path.exists(metadata_file):
            with open(metadata_file,'w+') as file:
                for main_id in bigg_metadata:
                    line = [main_id,'|']
                    line.append(f'enzyme_ec:{main_id}')
                    for db_id in bigg_metadata[main_id]:
                        line.append(f'bigg:{db_id}')
                    file.write('\t'.join(line)+'\n')
        os.remove(bigg2refs_file)

    def workflow_function(self):

        uniprot_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz'
        compressed_uniprot = f'{self.work_dir}{SPLITTER}uniprot_file.xml.gz'
        uncompressed_uniprot = f'{self.work_dir}{SPLITTER}uniprot_file.xml'
        hmm_file=f'{self.work_dir}{SPLITTER}uniprot_ec.hmm'

        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)

        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)

        if not os.path.exists(compressed_uniprot) or not os.path.exists(uncompressed_uniprot):
            self.download_file_ftp(uniprot_url, compressed_uniprot)

        if not os.path.exists(uncompressed_uniprot):
            self.unpack_gz(compressed_uniprot, uncompressed_uniprot)

        if not os.listdir(self.fasta_dir):
            dict_accession_seqs, dict_ec_accession = self.uniprot_xml_parser(uncompressed_uniprot)
            self.fasta_writer(dict_accession_seqs, dict_ec_accession)

        self.launch_fastas_msa()
        self.launch_aln_hmmer()
        self.merge_profiles(output_file=hmm_file)
        print(f'Finished generating {hmm_file}')


class HMM_Updater_Uniprot_Rhea(HMM_Updater):
    def __init__(self,work_dir,remove_files,min_seqs):
        HMM_Updater.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs)
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


        rhea2xrefs_file = f'{self.work_dir}{SPLITTER}rhea2xrefs.tsv'
        rhea2xrefs_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2xrefs.tsv'

        rhea2bigg_file=f'{self.work_dir}{SPLITTER}bigg_models_reactions.txt'
        rhea2bigg_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'


        if not os.path.exists(rhea2xrefs_file):
            self.download_file_ftp(rhea2xrefs_url, rhea2xrefs_file)

        if not os.path.exists(rhea2bigg_file):
            self.download_file_ftp(rhea2bigg_url, rhea2bigg_file)

        wanted_db='rhea'
        metadata_file = f'{self.work_dir}{SPLITTER}uniprot_{wanted_db}.tsv'
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
                if uniprot_id in seq_ids and len(seq_ids)>=self.min_seqs:
                    fasta_file = f'{self.fasta_dir}{main_id}.faa'
                    with open(fasta_file, 'a+') as file:
                        outline = f'>{uniprot_id}\n{sequence}\n'
                        file.write(outline)


    def workflow_function(self):

        uniprot_fastas_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
        compressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml.gz'
        uncompressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml'

        rhea2uiniprot_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv'
        rhea2uiniprot_file = f'{self.work_dir}{SPLITTER}rhea2uniprot.tsv'




        hmm_file=f'{self.work_dir}{SPLITTER}uniprot_rhea.hmm'


        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)

        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)


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


class HMM_Updater_Uniprot_Reactome(HMM_Updater):
    def __init__(self,work_dir,remove_files,min_seqs):
        HMM_Updater.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs)
        self.workflow_function()
        self.write_metadata()


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
                    if master_id not in res: res[master_id] = set()
                    res[master_id].add(uniprot_id)
                line = file.readline()
        return res



    def write_metadata(self):
        bigg2refs_file=f'{self.work_dir}{SPLITTER}bigg_models_reactions.txt'
        bigg2refs_file_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'

        if not os.path.exists(bigg2refs_file):
            self.download_file_ftp(bigg2refs_file_url, bigg2refs_file)

        wanted_db='reactome'
        metadata_file = f'{self.work_dir}{SPLITTER}uniprot_{wanted_db}.tsv'
        bigg_metadata=self.parse_bigg(bigg2refs_file,wanted_dbs=[wanted_db])
        if not os.path.exists(metadata_file):
            with open(metadata_file,'w+') as file:
                for main_id in bigg_metadata:
                    line = [main_id,'|']
                    line.append(f'{wanted_db}:{main_id}')
                    for db_id in bigg_metadata[main_id]:
                        line.append(f'bigg:{db_id}')
                    file.write('\t'.join(line)+'\n')
        os.remove(bigg2refs_file)

    def fasta_writer(self,uniprot_mapping, uncompressed_uniprot_fastas):
        uniprot_seqs=self.read_protein_fasta_generator(uncompressed_uniprot_fastas)
        for uniprot_seq in uniprot_seqs:
            uniprot_id,sequence=uniprot_seq
            for main_id in uniprot_mapping:
                seq_ids=uniprot_mapping[main_id]
                if uniprot_id in seq_ids and len(seq_ids)>=self.min_seqs:
                    fasta_file = f'{self.fasta_dir}{main_id}.faa'
                    with open(fasta_file, 'a+') as file:
                        outline = f'>{uniprot_id}\n{sequence}\n'
                        file.write(outline)


    def workflow_function(self):

        uniprot_fastas_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
        compressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml.gz'
        uncompressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml'

        reactome2uniprot_url='https://reactome.org/download/current/UniProt2ReactomeReactions.txt'
        reactome2uniprot_file = f'{self.work_dir}{SPLITTER}UniProt2ReactomeReactions.txt'

        hmm_file=f'{self.work_dir}{SPLITTER}uniprot_reactome.hmm'


        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)

        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)


        if not os.path.exists(compressed_uniprot_fastas) or not os.path.exists(uncompressed_uniprot_fastas):
            self.download_file_ftp(uniprot_fastas_url, compressed_uniprot_fastas)


        if not os.path.exists(reactome2uniprot_file):
            self.download_file_ftp(reactome2uniprot_url, reactome2uniprot_file)

        if not os.path.exists(uncompressed_uniprot_fastas):
            self.unpack_gz(compressed_uniprot_fastas, uncompressed_uniprot_fastas)

        if not os.listdir(self.fasta_dir):
            reactome_uniprot = self.parse_reactome2uniprot(reactome2uniprot_file)
            self.fasta_writer(reactome_uniprot, uncompressed_uniprot_fastas)


        self.launch_fastas_msa()
        self.launch_aln_hmmer()
        self.merge_profiles(output_file=hmm_file)
        print(f'Finished generating {hmm_file}')
        os.remove(reactome2uniprot_file)
        os.remove(compressed_uniprot_fastas)
        os.remove(uncompressed_uniprot_fastas)


class HMM_Updater_Uniprot_BIGG(HMM_Updater,Web_Connector):
    def __init__(self,work_dir,remove_files,min_seqs):
        HMM_Updater.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs)
        Web_Connector.__init__(self)
        self.mp_results = self.manager.list()
        self.workflow_function()
        self.write_metadata()


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
        print(f'Getting info for model {model_id} and gene {gene_id}')
        json_page=self.get_url_json(models_url)
        protein_sequence=json_page['protein_sequence']
        reactions=json_page['reactions']
        reactions_bigg=set()
        for i in reactions:
            bigg_id=i['bigg_id']
            reactions_bigg.add(bigg_id)
        return [gene_id,protein_sequence,reactions_bigg]

    def gene_info_worker_function(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            arg1, arg2 = record
            self.mp_results.append(self.get_gene_info(arg1, arg2))

    def launch_reaction_info_retrieval(self, model_id, genes_list):
        for gene_id in genes_list:
            self.queue.append([model_id, gene_id])
        self.processes_handler(self.gene_info_worker_function)
        while self.mp_results:
            yield self.mp_results.pop(0)

    def export_to_fasta(self,reaction_id,gene_id,protein_sequence):
        fasta_path = f'{self.fasta_dir}{SPLITTER}{reaction_id}.faa'
        with open(fasta_path, 'a+') as file:
            outline = f'>{gene_id}\n{protein_sequence}\n'
            file.write(outline)

    def fasta_writer(self):
        for model_id in self.get_all_models():
            print(f'Getting info for model {model_id}')
            genes_list=self.get_genes_model(model_id)
            reactions_generator=self.launch_reaction_info_retrieval(model_id,genes_list)
            for gene_id, protein_sequence, reactions_bigg in reactions_generator:
                for reaction_id in reactions_bigg:
                    self.export_to_fasta(reaction_id,gene_id,protein_sequence)
        for fasta_file in os.listdir(self.fasta_dir):
            fasta_path=f'{self.fasta_dir}{SPLITTER}{fasta_file}'
            len_fasta=self.get_seqs_count(fasta_path)
            if len_fasta<self.min_seqs:
                os.remove(fasta_path)

    def write_metadata(self):
        bigg_file=f'{self.work_dir}{SPLITTER}bigg_models_reactions.txt'
        bigg_url='http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
        if not os.path.exists(bigg_file):
            self.download_file_ftp(bigg_url, bigg_file)
        metadata_file = f'{self.work_dir}{SPLITTER}bigg.tsv'
        bigg_metadata=self.parse_bigg(bigg_file)
        reactions_ids=[i.replace('.hmm','') for i in os.listdir(self.hmm_dir)]
        if not os.path.exists(metadata_file):
            with open(metadata_file,'w+') as file:
                for main_id in reactions_ids:
                    if main_id in bigg_metadata:
                        line = [main_id,'|']
                        for db in bigg_metadata[main_id]:
                            for db_id in bigg_metadata[main_id]:
                                line.append(f'{db}:{db_id}')
                        file.write('\t'.join(line)+'\n')
        os.remove(bigg_file)




    def workflow_function(self):

        hmm_file=f'{self.work_dir}{SPLITTER}bigg.hmm'

        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)

        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)

        if not os.listdir(self.fasta_dir):
            self.fasta_writer()

        self.launch_fastas_msa()
        self.launch_aln_hmmer()
        self.merge_profiles(output_file=hmm_file)
        print(f'Finished generating {hmm_file}')


if __name__ == '__main__':
    print('Executing command:\n', ' '.join(argv))
    parser = argparse.ArgumentParser(description='An HMM generator tool using Uniprot sequences as reference and Rhea or EC to cluster these sequences\n'
                                     , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-db','--database', help='[required]\tClustering ID',choices=['ec', 'rhea','reactome','bigg'])
    parser.add_argument('-o', '--output_folder', help='[required]\tDirectory to save HMMs in')
    parser.add_argument('-ms', '--min_seqs',help='[optional]\tMinimum sequences per HMM. Default is 10')
    parser.add_argument('-rf', '--remove_files', action='store_true',help='[optional]\tuse this to remove files from previous runs.')

    args = parser.parse_args()
    database = args.database
    output_folder = args.output_folder
    min_seqs = args.min_seqs
    remove_files = args.remove_files
    if min_seqs:    min_seqs=int(min_seqs)
    else:           min_seqs=10
    if not output_folder:
        print('Missing output folder!')
    elif database=='rhea':
        updater = HMM_Updater_Uniprot_Rhea(output_folder, remove_files,min_seqs)
    elif database=='reactome':
        updater=HMM_Updater_Uniprot_Reactome(output_folder,remove_files,min_seqs)
    elif database=='ec':
        updater=HMM_Updater_Uniprot_EC(output_folder,remove_files,min_seqs)
    elif database=='bigg':
        updater=HMM_Updater_Uniprot_BIGG(output_folder,remove_files,min_seqs)
    else:
        print('Command is not valid')
