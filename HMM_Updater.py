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

__author__ = "Pedro Queirós, Polina Novikova"
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
        self.metadata_file = f'{self.work_dir}{SPLITTER}metadata.tsv'


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
            self.hmm_presser(hmm_file)

    def hmm_presser(self,hmm_file):
        software = 'hmmpress'
        command = f'{software} {hmm_file}'
        print('Running command:', command)
        subprocess.run(command, shell=True)


class HMM_Updater_Uniprot_EC(HMM_Updater):
    def __init__(self,work_dir,remove_files,min_seqs):
        HMM_Updater.__init__(self,work_dir=work_dir,remove_files=remove_files,min_seqs=min_seqs)
        self.workflow_function()


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

    def write_metadata(self,rhea2xrefs_path):
        if not os.path.exists(self.metadata_file):
            print('Parsing rhea2xrefs')
            rhea2ids={}
            with open(rhea2xrefs_path) as file:
                line=file.readline()
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    if line:
                        rhea_id,direction,master_id,db_id,db_type=line.split('\t')
                        if master_id not in rhea2ids: rhea2ids[master_id]={}
                        if db_type=='EC': db_type='enzyme_ec'
                        elif db_type=='METACYC': db_type='biocyc'
                        elif db_type=='ECOCYC': db_type='biocyc'
                        elif db_type=='KEGG_REACTION': db_type='kegg'
                        elif db_type=='GO':
                            db_type='go'
                            db_id=db_id.strip('GO:')

                        elif db_type=='REACTOME': db_type=None
                        elif db_type=='MACIE': db_type=None
                        if db_type:
                            if db_type not in rhea2ids[master_id]: rhea2ids[master_id][db_type]=set()
                            rhea2ids[master_id][db_type].add(db_id)
                    line=file.readline()
            with open(self.metadata_file,'w+') as file:
                for rhea_id in rhea2ids:
                    line = [rhea_id,'|']
                    for db_type in rhea2ids[rhea_id]:
                        for db_id in rhea2ids[rhea_id][db_type]:
                            line.append(f'{db_type}:{db_id}')
                    file.write('\t'.join(line)+'\n')

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

    def fasta_writer(self,rhea_uniprot, uncompressed_uniprot_fastas):
        uniprot_seqs=self.read_protein_fasta_generator(uncompressed_uniprot_fastas)
        for uniprot_seq in uniprot_seqs:
            uniprot_id,sequence=uniprot_seq
            for rhea_id in rhea_uniprot:
                seq_ids=rhea_uniprot[rhea_id]
                if uniprot_id in seq_ids and len(seq_ids)>=self.min_seqs:
                    fasta_file = f'{self.fasta_dir}{rhea_id}.faa'
                    with open(fasta_file, 'a+') as file:
                        outline = f'>{uniprot_id}\n{sequence}\n'
                        file.write(outline)

    def workflow_function(self):

        uniprot_fastas_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
        compressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml.gz'
        uncompressed_uniprot_fastas = f'{self.work_dir}{SPLITTER}uniprot_file.xml'

        rhea2uiniprot_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv'
        rhea2uiniprot_file = f'{self.work_dir}{SPLITTER}rhea2uniprot.tsv'

        rhea2xrefs_file = f'{self.work_dir}{SPLITTER}rhea2xrefs.tsv'
        rhea2xrefs_url='https://ftp.expasy.org/databases/rhea/tsv/rhea2xrefs.tsv'



        hmm_file=f'{self.work_dir}{SPLITTER}uniprot_rhea.hmm'


        if os.path.exists(self.work_dir) and self.remove_files:
            shutil.rmtree(self.work_dir)

        for directory in [self.work_dir, self.fasta_dir, self.aln_dir, self.hmm_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)


        if not os.path.exists(compressed_uniprot_fastas) or not os.path.exists(uncompressed_uniprot_fastas):
            self.download_file_ftp(uniprot_fastas_url, compressed_uniprot_fastas)


        if not os.path.exists(rhea2uiniprot_file):
            self.download_file_ftp(rhea2uiniprot_url, rhea2uiniprot_file)

        if not os.path.exists(rhea2xrefs_file):
            self.download_file_ftp(rhea2xrefs_url, rhea2xrefs_file)


        if not os.path.exists(uncompressed_uniprot_fastas):
            self.unpack_gz(compressed_uniprot_fastas, uncompressed_uniprot_fastas)

        if not os.listdir(self.fasta_dir):
            rhea_uniprot = self.parse_rhea2uniprot(rhea2uiniprot_file)
            self.fasta_writer(rhea_uniprot, uncompressed_uniprot_fastas)

        self.launch_fastas_msa()
        self.launch_aln_hmmer()
        self.merge_profiles(output_file=hmm_file)
        self.write_metadata(rhea2xrefs_file)
        print(f'Finished generating {hmm_file}')


if __name__ == '__main__':
    print('Executing command:\n', ' '.join(argv))
    parser = argparse.ArgumentParser(description='An HMM generator tool using Uniprot sequences as reference and Rhea or EC to cluster these sequences\n'
                                     , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-db','--database', help='[required]\tClustering ID',choices=['ec', 'rhea'])
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
    elif database=='ec':
        updater=HMM_Updater_Uniprot_EC(output_folder,remove_files,min_seqs)
    else:
        print('Command is not valid')
