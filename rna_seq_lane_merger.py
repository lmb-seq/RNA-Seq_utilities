import os, glob, argparse, re, subprocess, pandas
from cell_bio_util import cell_bio_util as util


def glob_lister(submission_form):
  '''Reads in the CRUKCI_SLX_Submission.xlsx file and returns a list of the
  glob variables to be used by the "globber" function (The "Index" column) in
  the file.
  
  Parameters
  ----------
  filepath (string / os.path):
    File location for glob variables. Preferably full path.
  
  Returns
  -------
  index_list (List):
    Containing glob variables.
  
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
  excel_file_find_start = excel_file[excel_file[0] == 'Name']
  excel_file_start = int(excel_file_find_start.index[0]) + 1

  samples_information = pandas.read_excel(submission_form, header = excel_file_start, usecols = ['Index'])
  samples_information['Index'] = samples_information['Index'].str.replace('-', '_')

  index_list = list(samples_information['Index'])

  return(index_list)


def globber(directory, glob_list_pattern):
  '''Takes in the directory to search within and a list of glob variables (e.g.
  from the "glob_lister" function). Returns a dictionary where the keys are the
  glob patterns and the values are the files associated with the glob search
  within the specified directory.

  Parameters
  ----------
  directory (string / os.path):
    Full path to perform glob search in. 
  pair_tags (list):
    Same as PRAGUI's pair tags input.
  
  Returns
  -------
  Dictionary;
    Keys (string):
      The glob patterns input ID (pair tag if applicable).
    Values: (list):
      Files to be merged.  
  
   
  '''

  index_files_dict = {}

  util.info('Obtaining file paths for for each of the listed indexes')
  for i in glob_list_pattern:
    glob_pattern = os.path.join(directory, '*{0}*.fq.gz'.format(i))
    globbing = glob.glob(glob_pattern)
    index_files_dict[i] = globbing
    util.info('Obtained file paths for files associated with {0}'.format(i))

  return(index_files_dict)


def lane_merger_preparation(indexed_files, paired_single):
  '''If the input files for PRAGUI are paired reads then this function will
  allow these pairs to be treated separately. e.g. files "FileA-r_1" and
  "FileA-r_2.fq.gz" will be merged separately. 
  
  Parameters
  ----------
  indexed_files (dictionary)
    Output of "globber" function. 
  pair_tags (list)
    Same as PRAGUI's pair tags input.
  
  Returns
  -------
  Dictionary:
    Keys (string):
      The glob patterns input (same as with the "globber" function")
    Values: (list):
      Files to be merged.
  
  
  '''

  util.info('Beginning lane merger preparation')
  laned_files = {}
  for indexes in indexed_files:
    if paired_single == 'paired':
      util.info('Separating files for {0} based on pair tags'.format(indexes))
      for tags in paired_tags:
        re_pattern_pairs = re.compile('.*?{0}.fq.gz'.format(tags))
        files = indexed_files[indexes]
        pair_seperated_list = list(filter(re_pattern_pairs.match, files))
        laned_files['{0} {1}'.format(indexes, tags)] = pair_seperated_list
        util.info('{0} files found for index {1}, read pair {2}'.format(len(pair_seperated_list), indexes, tags))
    elif paired_single == 'single':
      util.info('No pair tags defined, including all files associated with {0}'.format(indexes))
      laned_files[indexes] = indexed_files[indexes]
      util.info('{0} files found for index {1}'.format(len(indexed_files[indexes]), indexes))

  return(laned_files)


def lane_merged_subfolder(working_directory):
  '''Creates a "lane_merged" subfolder within the working directory. 
  
  Parameters
  ----------
  working_directory (string / os.path):
    Location of input files.
  
  Returns
  -------
  Boolean (True / False):
    Can user place files within created subfolder?
  subfolder (string / os.path):
    Full path of created subfolder.
  
  '''

  subfolder = os.path.join(working_directory, 'lane_merged')
  util.info('Attempting to create "lane_merged_files" folder within {0}'.format(working_directory))
  try:
    os.mkdir(subfolder)
  except(OSError):
    if os.path.isdir(subfolder):
      util.info('Directory already exists: {0}'.format(subfolder))
      return(True, subfolder)
    else:
      util.critical('Unable to create directory: {0}. Cannot proceed further'.format(subfolder))
      return(False, subfolder)
  else:
    util.info('Subfolder creation successful: {0}'.format(subfolder))
    return(True, subfolder)


def merged_filename(file_list, lane_tags, output_subfolder):
  '''Takes the list of files to merge and provides an output file (filename).
  
  Parameters
  ----------
  file_list (List):
    Full paths of input files.
  lane_tags (List):
    Tags that identify samples' lane e.g. "s_1 s_2".
  output_subfolder (string / os.path):
    Location for output file.
  
  
  Returns
  -------
  output_file_str (string / os.path):
    Full path of output file to be created by this script.
  
  '''

  output_files = []
  for tags in lane_tags:
    for full_file_path in file_list:
      if(tags in full_file_path):

        path_split = os.path.split(full_file_path)

        new_file = path_split[1].replace(tags, 'merged')
        new_full_file_path = os.path.join(output_subfolder, new_file)
        output_files.append(new_full_file_path)

  output_file_str = list(set(output_files))[0]
  return(output_file_str)


def lane_merger(working_directory, files_to_merge, lane_tags):
  '''Performs the merging of the input files. zcat to read in, pigz to create
  the merged file - multi-threaded alternative to gzip for faster compression. 
  
  Parameters
  ----------
  working_directory (string / os.path)
    Location of input files.
  files_to_merge (Dictionary):
    Keys (string):
      The glob patterns input (same as with the "globber" function")
    Values: (list):
      Files to be merged.
  lane_tags (List):
    Tags that identify samples' lane e.g. "s_1 s_2".
  
  
  Returns
  -------
  No return
  
  '''

  subfolder_check, subfolder = lane_merged_subfolder(working_directory)

  if subfolder_check == False:
    util.critical('Terminating script early')

  util.info('Beginning lane merger for files')

  for index_files in files_to_merge:
    util.info('Merging {0}'.format(index_files))
    input_files = files_to_merge[index_files]
    output_file_name_pre = merged_filename(input_files, lane_tags, subfolder)
    output_file_name = os.path.join(subfolder, output_file_name_pre)
    with open(output_file_name, 'wb') as outfile:
      zcat_files = util.run(['zcat'] + input_files, stdout = subprocess.PIPE)
      pigz_output = util.run(['pigz', '-c'], stdin = zcat_files.stdout, stdout = outfile)
      zcat_files.wait()
    util.info('Output file {0} created'.format(output_file_name))
  util.info('All lane files merged')


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description = 'Merge RNA-Seq files across lanes.')
  parser.add_argument('-f', '--submission_form',
                      type = str,
                      required = True,
                      metavar = 'FILENAME',
                      help = 'Path to the submission form provided (e.g. CRUKCI_SLX_Submission.xlsx) - Please provide full path and ensure this file is in same folder as the RNA-Seq files.')

  parser.add_argument('-l', '--lane_tags',
                      type = str,
                      required = True,
                      nargs = '+',
                      help = 'Tags that identify samples\' RNA-Seq lanes e.g. "s_1 s_2".')

  group = parser.add_mutually_exclusive_group(required = True) # Sets --single_end and --paired_end as mutually exclusive arguments
  group.add_argument('-s', '--single_end', action = 'store_true',
                     help = 'Flag if RNA-Seq data are single end reads. Mutually exclusive with the -p/--paired_end argument.')
  group.add_argument('-p', '--paired_end', nargs = 2, metavar = '<PAIR_TAG>',
                     help = 'Flag if RNA-Seq data are paired end reads. Mutually exclusive with the -s/--single_end argument. Provide pair tags, this will be the same as PRAGUI\'s "pair_tags" argument.')

  args = parser.parse_args()
  if args.single_end == True:
    paired_single = 'single'
    paired_tags = None
  else:
    paired_single = 'paired'
    paired_tags = args.paired_end

  working_directory = os.path.split(args.submission_form)[0]

  glob_list = glob_lister(args.submission_form)
  indexed_files = globber(working_directory, glob_list)
  files_to_merge = lane_merger_preparation(indexed_files, paired_single)
  merged_files = lane_merger(working_directory, files_to_merge, args.lane_tags)
  util.info('Process complete')
