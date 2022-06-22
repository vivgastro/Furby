#!/home/observer/miniconda2/bin/python

import numpy as np
import os, sys, warnings, argparse

from sigpyproc.Readers import FilReader as F
from Furby_reader import Furby_reader as Fr

def get_injected_data(furby, filt, samp):
  ff = Fr(furby)
  ff_data = ff.read_data()
  filt_data = filt.readBlock(samp, ff.header.NSAMPS)
  rms_of_filt_data_per_chan = filt_data.std(axis=1)
  added = filt_data + (ff_data * rms_of_filt_data_per_chan[:, None])
  return added.astype(filt.header.dtype)

def write_to_filt(data, out):
  if data.dtype != out.dtype:
     warnings.warn("Given data (dtype={0}) will be unasfely cast to the requested dtype={1} before being written out".format(data.dtype, o.dtype), stacklevel=1)
  out.cwrite(data.T.ravel().astype(out.dtype, casting='unsafe'))

def copy_filts(inp, out, start, end, gulp=8192):
  for nsamps, ii, d in inp.readPlan(gulp, start=start, nsamps=end-start, verbose=False):
    write_to_filt(d, out)

def assert_isamp_sep(isamps):
  x = [0]
  x.extend(list(isamps))
  x = np.array(x)
  diff = x[1:] - x[:-1]
  if np.any(diff < 9000):
    raise ValueError("Injection time stamps cannot be less than 9000 samples apart")

def main():
  filt = args.filt
  furbies = args.furbies
  isamps = np.array(args.isamps)
  assert len(isamps) == len(furbies)

  f = F(filt)
  isamps = sorted((isamps/f.header.tsamp).astype(np.int))
  assert_isamp_sep(isamps)
  o = f.header.prepOutfile("injected_{}".format(f.filename))
  for ii, isamp in enumerate(isamps):
    if ii ==0:
      copy_filts(inp = f, out = o, start= 0, end = isamp)
    
    injected_data = get_injected_data(furby = furbies[ii], filt = f, samp = isamp)
    write_to_filt(data = injected_data, out=o)
    sys.stdout.write("Injected {0} in {1} at {2:.2f}\n".format(furbies[ii], f.filename, isamp*f.header.tsamp))
    
    if ii == len(isamps)-1:
      copy_filts(inp=f, out=o, start = isamp+Fr(furbies[ii]).header.NSAMPS, end = f.header.nsamples)
    else:
      copy_filts(inp=f, out=o, start = isamp+Fr(furbies[ii]).header.NSAMPS, end = isamps[ii+1])

if __name__ == '__main__':
  a = argparse.ArgumentParser()
  a.add_argument('-filt', type=str, help="Filterbank to inject in")
  a.add_argument('-furbies', type=str, nargs='+', help="Path to furby files")
  a.add_argument('-isamps', type=float, nargs='+', help="Injection time stamps in seconds")
  
  args = a.parse_args()
  main()
