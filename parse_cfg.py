import sys, os
from check_type import check_type as ctp
from collections import namedtuple

def parse_cfg(cfg_file, return_dict = False, comment_identifier = ("#"), delim = None, ret_comments=False, ret_studs=False):
  '''
  Reads any text configuration file into a namedtuple (or dict).

  The first column becomes the key and the second/rest column(s) become(s) the value. If there are more than one values for a key, all the values stay in a string as a single value, unless 'delim' is specified in which case the values will be split into an array using 'delim' as the delimeter.
  (Note: Any occurence of the specified 'delim' in the first column (i.e. keys) will be ignored)

  You may specify a comment_identifier string, which by default is '#'. All lines starting with the 'comment_identifier' are ignored unless ret_comments is true when all the comments are appended to a list with key name 'Comments_'. Empty comments will be ignored. Any line which has a comment towards the end (e.g. :  x  5  #sample comment) will be ignored as well and the line will be processed as if it had no comment.

  Any lines with only one column will be ignored unless ret_studs is True when they are appended to a list with key name 'Studs_'.

  Only alphanumeric characters are supported as valid first characters of keys if returning a named_tuple (default). Any string will work as key if return_dict is True.
  '''
  if not os.path.exists(cfg_file):
    raise IOError("File {0} does not exist".format(cfg_file))
  
  if not os.access(cfg_file, os.R_OK):
    raise IOError("Do not have read permissions on {0}".format(cfg_file))
  
  conf_dict = {}
  if ret_comments:
    conf_dict["Comments_"] = []
  if ret_studs:
    conf_dict["Studs_"] = []

  with open(cfg_file) as c:
    while True:
      line = c.readline()
      if line=="":
        break
      elif line.strip() == "":
        continue
      elif len(line.split()) < 2:
        if ret_studs:
          conf_dict["Studs_"].append(line.strip())
        continue
      elif line.strip().startswith(comment_identifier):
        if len(line.strip().strip(delim)) > 0 and ret_comments:
          conf_dict["Comments_"].append(line.strip())
        continue
      else:
        key = line.split()[0].strip()
        val = line.strip()[len(key):].strip().split(comment_identifier)[0].strip()
        val = ctp(val)
        if delim is not None and isinstance(val, str):
          x = val.split(delim)
          if len(x) > 1:
            x = [ctp(i) for i in x]
            val = x
        conf_dict[key] = val
  
  if return_dict:
    return conf_dict.copy()

  tmp = namedtuple("CONF", conf_dict.keys())
  parsed_namedtuple = tmp(*conf_dict.values())

  return parsed_namedtuple

