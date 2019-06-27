import os, glob, argparse, sys, pandas, ftplib, hashlib

module_path = os.path.realpath(__file__)
rnaseq_utilities_directory = os.path.dirname(module_path)
utilities_directory = os.path.split(rnaseq_utilities_directory)[0]
sys.path.append(utilities_directory)

from cell_bio_util import cell_bio_util as util


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

  directory_absolute_path = os.path.abspath(directory)

  if os.path.isdir(directory_absolute_path):
    return(directory_absolute_path)
  else:
    message = 'Invalid --{0} location. Please ensure directory exists before proceeding:\n\t{1}'.format(check_type, directory)
    util.critical(message)


def glob_lister(submission_form):
  '''Reads in the CRUKCI_SLX_Submission.xlsx file and returns the SLX number and fq.gz file prefixes.
  
  Parameters
  ----------
  filepath (string / os.path):
    File location for glob variables. Preferably full path.
  
  Returns
  -------
  samples_information (pandas dataframe):
    Pandas dataframe containing the file index (prefix)
  
  slx_id (string):
    SLX ID from CRUKCI
  
  '''

  file_check = os.path.isfile(submission_form)
  if file_check == True:
    util.info('Submission form found at {0}'.format(submission_form))
  else:
    error_message = 'Submission form not found, please ensure this file is in\n-->\t{0}'.format(working_directory)
    util.critical(error_message)

  util.info('Reading in index file')
  excel_file = pandas.read_excel(submission_form) # Find the start row to skip passed the header information
  util.info('Index file read')
  slx_id_raw = excel_file[excel_file[1] == 'SLX Identifier']
  slx_row = slx_id_raw.index[0]
  slx_id = 'SLX-{0}'.format(excel_file.iloc[slx_row]['Unnamed: 2'])

  excel_file_find_start = excel_file[excel_file[1] == 'Name']
  excel_file_start = int(excel_file_find_start.index[0]) + 1

  samples_information = pandas.read_excel(submission_form, header = excel_file_start, usecols = ['Name', 'Index'])
  samples_information['Index'] = samples_information['Index'].str.replace('-', '_')

  return(samples_information, slx_id)


def ftp_server_connection():
  ''' Establishes connection to the CRUK FTP server.
  This version uses the de Bono lab's user log credentials
  
  Returns
  -------
  ftp_server (ftblib object):
  
  '''

  ftp_server_url = 'ftp1.cruk.cam.ac.uk'
  util.info('Accessing FTP server {0}'.format(ftp_server_url))
  ftp_server = ftplib.FTP(ftp_server_url)
  ftp_server.login(user = 'lmb_debono', passwd = 'shiphill50')
  util.info('Logged into FTP server')

  return(ftp_server)


def ftp_download_files(ftp_server, slx_id):
  '''Downloads the fastq (.fq.gz files).
  
  Parameters
  ----------
  ftp_server (ftblib object):
    connection to FTP server
  
  slx_id (string)
    The user's SLX ID for the wanted files
  
  Returns
  -------
  downloaded_files (list):
    List of downloaded files (absolute path)_
  
  
  '''

  file_prefix = slx_id
  ftp_list_dir = ftp_server.nlst()
  downloaded_files = []

  util.info('Downloading files beginning with {0}'.format(file_prefix))

  for file in ftp_list_dir:
    if file.startswith(file_prefix):
      file_path = os.path.join(working_directory, '{0}'.format(file))
      util.info('Attempting to download file {0}'.format(file))

      if os.path.isfile(file_path) == True:
        util.info('File already exists, skipping'.format(file))
      else:
        with open(file_path, 'wb') as download_file:
          ftp_server.retrbinary('RETR ' + file, download_file.write, 1024)
        util.info('File downloaded to {0}'.format(file_path))

      downloaded_files.append(file_path)

  return(downloaded_files)


def file_md5_check(downloaded_files_list):
  '''Performs a check to see if MD5 checksum matches with expected MD5 hash
  
  Parameters
  ----------
  downloaded_files_list (list):
    List containing files that were downloaded
  
  Returns
  -------
  Boolean (True/False)
    Does expected MD5 hash values (within CRUK provided .md5sums.txt files match with computed MD5 checksums?
  
  '''

  list_directory = downloaded_files_list
  md5_hash_value_files = []
  md5_check_hash_dictionary = {}
  failed_downloads = []

  for file in list_directory:
    if file.endswith('.md5sums.txt'):
      md5_hash_value_files.append(file)

  for md5_checksum_file in md5_hash_value_files:
    md5_checksum_file_path = os.path.join(working_directory, md5_checksum_file)
    with open(md5_checksum_file_path, 'r') as check_file:
      file_read = check_file.readlines()
      for line in file_read:
        md5_hash, file = line.strip().split('  ')
        md5_check_hash_dictionary[file] = md5_hash

  for sample_file in md5_check_hash_dictionary:

    file_path = os.path.join(working_directory, sample_file)
    file_hashing = hashlib.md5()
    with open(file_path, 'rb') as file_to_check:
      for chunks in iter(lambda: file_to_check.read(4096), b""):
        file_hashing.update(chunks)

    util.info('Performing MD5 hash check')
    hashing_value_calculated = file_hashing.hexdigest()
    hashing_value_expected = md5_check_hash_dictionary[sample_file]
    util.info('Computed MD5 hash value = {0}'.format(hashing_value_calculated))
    util.info('Expected MD5 hash value = {0}'.format(hashing_value_expected))

    if hashing_value_calculated == hashing_value_expected:
      util.info('Computed and expected MD5 hash values match')
    else:
      util.warn('Computed and expected MD5 hash values do not match')
      failed_downloads.append(sample_file)

  if len(failed_downloads) == 0:
    util.info('All files downloaded successfully')
    return(True)
  else:
    util.warn('{0} files did not pass MD5 checksum test'.format(len(failed_downloads)))
    for failed in failed_downloads:
      util.info('Deleting file {0}'.format(failed))
      failed_file = os.path.join(working_directory, failed)
      os.remove(failed_file)
    return(False)


def samples_csv_writer(working_directory, slx_id, samples_information):
  '''Automatically creates a samples.csv file based on what was downloaded and included in CRUKCI_SLX_Submission.xlsx file
  
  Parameters
  ----------
  working_directory (os.path / string)
    Allows samples.csv to be filled with fq.gz files' absolute path
  
  slx_id
    Allows relevant files to be IDed using file prefix
  
  samples_information
    Allows relevant and specific files to be IDed using file midfix
    
  '''

  output_file = os.path.join(working_directory, 'samples.csv')
  util.info('Writing to {0}'.format(output_file))

  samples_information['read1'] = ''
  samples_information['read2'] = ''
  samples_information['condition'] = ''

  for index, row in samples_information.iterrows():
    glob_pattern = os.path.join(working_directory, '{0}.{1}*.r_*.fq.gz'.format(slx_id, row['Index']))
    globbing = sorted(glob.glob(glob_pattern))
    samples_information['read1'][index] = globbing[0]

    if len(globbing) == 2:
      samples_information['read2'][index] = globbing[1]

  del samples_information['Index']
  samples_information = samples_information.rename(columns = {'Name': 'samples'})

  samples_information.to_csv(output_file, sep = '\t', index = False)
  util.info('Write complete. See file located {0}'.format(output_file))


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description = 'Download RNA-Seq files from CRUK\'s FTP server')
  parser.add_argument('-f', '--submission_form',
                      type = str,
                      required = True,
                      metavar = 'FILENAME',
                      help = 'Path to the submission form (e.g. CRUKCI_SLX_Submission.xlsx) - Please provide full path and ensure this file is in same folder as where you wish to download the RNA-Seq files to.')

  args = parser.parse_args()

  working_directory = os.path.abspath(os.path.split(args.submission_form)[0])
  directory_check = check_directory(working_directory, 'working_directory')

  samples_information, slx_id = glob_lister(os.path.abspath(args.submission_form))
  ftp_server = ftp_server_connection()
  downloaded_files = ftp_download_files(ftp_server, slx_id)
  ftp_server.quit()

  retries = 0
  md5_check = False
  while md5_check == False:
    if retries >= 3:
      util.critical('{0} retries at downloading files have failed. Please try again later'.format(retries))

    md5_check = file_md5_check(downloaded_files)
    retries += 1

  samples_csv_writer(working_directory, slx_id, samples_information)

  util.info('\nrun complete')
