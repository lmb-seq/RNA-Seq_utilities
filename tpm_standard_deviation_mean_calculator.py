import argparse, os, sys, pandas, re, unicodedata
from cell_bio_util import cell_bio_util as util


def check_file(file_path):
  '''Checks to see if a passed file exits
  
  Parameters
  ----------
  file_path (string / os.path):
    File to be checked
    
  Returns
  -------
  file_path (string / os.path):
    Path to file returned if valid or exception if invalid
  
  '''

  file_check = os.path.isfile(file_path)
  if file_check == True:
    util.info('TPM file found')
    return(file_path)
  else:
    sys.exit('Unable to locate file check directory: {0}'.format(file_path))


def working_directory_finder(samples_tpm_path):
  '''Uses the TPM file and sets the working directory
  
  Parameters
  ----------
  samples_tpm_path (string / os.path):
    Location of TPM file
    
  Returns
  -------
  working_directory (string / os.path):
    Directory to the TPM file
    
  tpm_file_name (string / os.path):
    TPM file name
  
  '''

  working_directory, tpm_file_name = os.path.split(samples_tpm_path)
  return(working_directory, tpm_file_name)


def read_in_tpm(tpm_file_location):
  '''Reads in the information within TPM file using the pandas module
  
  Parameters
  ----------
  tpm_file_location (string / os.path):
    Location of TPM file
    
  Returns
  -------
  tpm_file (pandas csv object):
    Pandas object containing TMP data
  
  '''

  tpm_file = pandas.read_csv(tpm_file_location, delimiter = '\t')
  util.info('TPM file read in')
  return(tpm_file)


def gene_name_converter():
  '''Provides a reference to allow gene IDs to be converted human readable gene name abbreviations 
  (e.g. "WBGene00000001" to "aap-1" etc)
  This version of the script only deals with C. elegans
  
  Parameters
  ----------
  NONE
    
  Returns
  -------
  gene_id (pandas csv object):
    Pandas object containing wormbase geneIDs and corresponding gene name abbreviations
  
  '''

  gene_id_file = os.path.join(os.sep, 'data1', 'geneIDs', 'c_elegans.canonical_bioproject.current.geneIDs.txt')
  util.info('Reading in geneIDs from {0}'.format(gene_id_file))
  gene_id_csv = pandas.read_csv(gene_id_file, delimiter = ',', names = [0, 'geneName', 'gene', 'transcript_id', 4])
  gene_id = gene_id_csv.filter(['geneName', 'gene', 'transcript_id'])
  util.info('Gene names/IDs read in')

  gene_id = gene_id.set_index('geneName')
  util.info('Filling in missing gene names with transcript IDs if available')
  gene_id.loc[gene_id['gene'].isnull(), 'gene'] = gene_id['transcript_id']

  return(gene_id)


def slugify(value):
  '''Normalises string; converts to lowercase, removes non-alphanumeric characters and converts spaces to hyphens.
  This function was developed to remove the need to install external slugify libraries
  
  Parameters
  ----------
  value (string):
    string to be sligified
    
  Returns
  -------
  value (string):
    string that has been sligified
  
  '''

  value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
  value = (re.sub('[^\w\s-]', '', value.decode('UTF8')).strip().lower())
  value = re.sub('[-\s]+', '-', value)
  return(value)


def samples_file_conditions_finder(tpm_file):
  '''Finds the user set conditions for the RNA-Seq data
  Reads in original CSV file to obtain information
  
  Parameters
  ----------
  tpm_file (string / os.path):
    TPM file location
  
  Returns
  -------
  sample_conditions (dictionary)
    Dictionary object; Keys are (slugified) conditions and values are lists containing the sample names 
  
  '''

  sample_conditions = {}

  util.info('Reading in samples CSV file to obtain sample names and conditions')
  csv_file = os.path.join(working_directory, tpm_file.replace('_tpm.txt', ''))
  csv_read = pandas.read_csv(csv_file, delimiter = '\t')

  for row in csv_read.iterrows():
    sample, condition = (row[1]['samples'], row[1]['condition'])
    slugged_condition = slugify(condition)
    if not slugged_condition in sample_conditions:
      sample_conditions[slugged_condition] = []

    sample_conditions[slugged_condition].append(sample)

  util.info('Sample names and conditions read in')

  return(sample_conditions)


def output_file_creator(tpm_read, sample_conditions, gene_ids):
  '''Create output files (one for each condition)
  Each file contains the gene ID, gene abbreviation, the TPM values for each replicate and the TPM mean and standard deviations across the replicates
  
  Parameters
  ----------
  tpm_read (pandas object):
    Pandas object containing TPM data
  
  sample_conditions (dictionary):
    Dictionary object; Keys are (slugified) conditions and values are lists containing the sample names 
  
  gene_ids (pandas csv object):
    Pandas object containing wormbase geneIDs and corresponding gene name abbreviations
  
  Returns
  -------
  Output files (Microsoft Excel files)
    Files containing TPM values with calculated mean and standard deviations
    
  '''

  merged_data = pandas.merge(tpm_read, gene_ids, left_on = 'geneName', right_index = True)
  util.info('Merging data with gene name/ID references')

  number_of_files = len(sample_conditions)
  output_file_list = '\n'.join(sample_conditions)
  util.info('Creating {0} files; for each condition:\n{1}'.format(number_of_files, output_file_list))

  for condition in sample_conditions:
    conditioned_samples = sample_conditions[condition]
    new_dataframe = merged_data.filter(items = ['geneName', 'gene'] + conditioned_samples)

    util.info('Calculating mean average for {0}'.format(condition))
    new_dataframe['tpm_mean'] = new_dataframe.reindex(columns = conditioned_samples).mean(axis = 1)

    util.info('Calculating standard_deviation for {0}'.format(condition))
    new_dataframe['tpm_standard_deviation'] = new_dataframe.reindex(columns = conditioned_samples).std(axis = 1)

    new_file_name = 'TPM_std_dev_{0}.xlsx'.format(condition)
    new_file_path = os.path.join(working_directory, new_file_name)
    util.info('Creating output file for condition {0} as {1}'.format(condition, new_file_path))
    new_dataframe.to_excel(new_file_name)
    util.info('File saved')


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description = 'Generate mean and standard deviation for TPM values')
  parser.add_argument('-t', '--tpm_file',
                        type = str,
                        required = True,
                        metavar = 'FILENAME',
                        help = 'Full path for the samples.csv file of interest')

  user_args = parser.parse_args()
  tpm_file_path = check_file(user_args.tpm_file)

  working_directory, tpm_file_name = working_directory_finder(tpm_file_path)
  os.chdir(working_directory)

  read_tpm = read_in_tpm(tpm_file_path)
  gene_ids = gene_name_converter()
  sample_conditions = samples_file_conditions_finder(tpm_file_name)
  output_file_creator(read_tpm, sample_conditions, gene_ids)

  util.info('Task complete')
