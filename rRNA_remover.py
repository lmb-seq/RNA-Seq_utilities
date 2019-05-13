import os, subprocess, argparse, sys

current_path = os.path.realpath(__file__)
current_path = os.path.join(current_path, 'cell_bio_util')

sys.path.append(current_path)
import cell_bio_util as util


def check_directory(directory, check_type):
  '''Checks to see if a passed directory path exits
  
  Parameters
  ----------
  Directory (string / os.path):
    Directory to be checked
  
  Check_type (string)
    One of the argparse flags
  
  Returns
  -------
  Directory (string / os.path):
    Path returned if valid or exception if invalid
  
  '''

  if os.path.isdir(directory):
    return(directory)
  else:
    raise Exception('Invalid --{0} location. Please ensure directory exists before proceeding:\n\t{1}'.format(check_type, directory))


def check_rRNA_library(rRNA_genome_path):
  '''Checks to see if a passed rRNA library exits
  
  Parameters
  ----------
  rRNA_genome_path (string / os.path):
    Filepath to be checked
  
  Returns
  -------
  Directory (string / os.path):
    Path returned if valid or error if invalid - terminating this script
  
  '''

  directory = os.path.split(rRNA_genome_path)

  directory_checked = check_directory(directory[0], 'rRNA_library')
  list_files = os.listdir(directory_checked)
  list_checked_list = []
  for files in list_files:
    if files.startswith(directory[1]) and files.endswith('.bt2'):
      list_checked_list.append(True)
    else:
      list_checked_list.append(False)

  list_checked = set(list_checked_list)
  if True in list_checked:
    util.info('rRNA library path appears valid')
    return(rRNA_genome_path)
  else:
    util.error('Unable to locate rRNA library (bt2 files) within specified rRNA_library path')
    sys.exit()


def gzip_file_list(working_directory):
  '''Gets list of all .fq.gz files in the working directory (except .lostreads.fq.gz files)
  
  Parameters
  ----------
  working_directory (string / os.path):
    Directory to search for fq.gz files
  
  Returns
  -------
  fastq_gz_files(list)
    List object containing .fq.gz filenames (i.e. sample files)
  
  '''

  list_files = os.listdir(working_directory)
  fastq_gz_files = []

  for files in list_files:
    if files.endswith('.lostreads.fq.gz'):
      pass
    elif files.endswith('.fq.gz'):
      fastq_gz_files.append(files)
    else:
      continue

  if len(list_files) == 0:
    util.error('There are no gzipped fastq (FILENMAE.fq.gz) files within specified directory')
    sys.exit()
  util.info('List of gzipped fastq files read into script')

  return(fastq_gz_files)


def paired_reads_finder(fastq_gz_files):
  '''Sorts .fq.gz files so paired reads are grouped for later co-processing
  This version assumes CRUKCI file naming format
  
  Parameters
  ----------
  fastq_gz_files (list):
    List object containing .fq.gz filenames (i.e. sample files)
  
  Returns
  -------
  paired_reads (dictionary):
    Dictionary object; Keys are filename prefixes and values are the paired read files   
  
  '''

  paired_reads = {}
  for files in fastq_gz_files:
    root_file_name_split = files.split('.') # Assumes CRUKCI file naming format
    root_file_name = '.'.join(root_file_name_split[0:4])
    pair_number = root_file_name_split[4]

    if not root_file_name in paired_reads:
      paired_reads[root_file_name] = {'1': '', '2': ''}

    if pair_number == 'r_1':
      paired_reads[root_file_name]['1'] = files
    elif pair_number == 'r_2':
      paired_reads[root_file_name]['2'] = files

  return(paired_reads)


def output_preperation(working_directory):
  '''Creates folders for the output and log files
  
  Parameters
  ----------
  working_directory (string / os.path):
    Directory to place sub-folders
  
  Returns
  -------
  output_subdirectory (string / os.path):
    Path to the output files sub-folder 
 
  '''
  subfolder = 'rRNA_processed'
  output_subdirectory = os.path.join(working_directory, subfolder)

  if not os.path.exists(output_subdirectory):
    util.info('Creating sub-folder {0} within {1}'.format(subfolder, working_directory))
    os.mkdir(os.path.join(output_subdirectory))
  util.info('Newly generated files will be stored within {0}/\n'.format(os.path.join(output_subdirectory)))

  if not os.path.exists(os.path.join(output_subdirectory, 'log_files')):
    util.info('Creating log_files sub-folder within {0}'.format(output_subdirectory))
    os.mkdir(os.path.join(output_subdirectory, 'log_files'))

  return(output_subdirectory)


def rrna_removal(paired_reads, output_subdirectory):
  '''Runs the bowtie2 commands on the reads to remove rRNA data.
  Currently only tested on paired data, need to update for single end data.
  
  Parameters
  ----------
  paired_reads (dictionary):
    Dictionary object; Keys are filename prefixes and values are the paired read files

  output_subdirectory (string / os.path):
    Path to the output files sub-folder 
    
  '''

  pair_number = 0
  for entries in paired_reads:
    pair_number += 1
    util.info('Processing pair number {0} of {1}: {2}'.format(pair_number, len(paired_reads), entries))

    if os.path.exists(os.path.join(output_subdirectory, 'logs_{0}.txt'.format(entries))):
      util.warning('{0} already processed, skipping'.format(entries))
      continue

    command = ['bowtie2', '--phred33', '-D', '20', '-R', '3', '-N', '1', '-L', '20',
               '-i', 'S,1,0.50', '-x', args.rRNA_library,
               '-X', '1000', '--dovetail', '-1', paired_reads[entries]['1'], '-2', paired_reads[entries]['2'],
               '-S', os.path.join(output_subdirectory, 'ribo_aligns_{0}.sam'.format(entries)),
               '--un-conc-gz', os.path.join(output_subdirectory, '{0}_rRNA_processed_r_%.fq.gz'.format(entries)),
               '--np', '0']

    with open(os.path.join(os.sep, output_subdirectory, 'log_files', 'logs_{0}.txt'.format(entries)), 'w') as stdout_file:
      stdout_file.write('\n\nRead pair file prefix: {0}'.format(entries))
      subprocess.run(command, stdout = stdout_file, stderr = stdout_file)

    os.remove(os.path.join(output_subdirectory, 'ribo_aligns_{0}.sam'.format(entries))) # Need to delete these .sam files otherwise accumulation of many large files


parser = argparse.ArgumentParser(description = 'Remove rRNA reads for RNA-Seq data')
parser.add_argument('-d', '--directory', help = 'Specify the location of the RNA-Seq data', type = str, metavar = 'DIRECTORY', required = True)
parser.add_argument('-l', '--rRNA_library', help = 'Specify location of the rRNA genome library. Use bowtie2 "-x" style argument. Default is C. elegans library.',
                    type = str, metavar = 'FILE', default = '/home/paulafp/c_elegans_rDNA/c_elegans_concat_rDNA')

args = parser.parse_args()

working_directory = check_directory(args.directory, 'working_directory')
os.chdir(working_directory)

check_rRNA_library(args.rRNA_library)
fastq_gz_files = gzip_file_list(working_directory)
paired_reads = paired_reads_finder(fastq_gz_files)
output_subdirectory = output_preperation(working_directory)
rrna_removal(paired_reads, output_subdirectory)
