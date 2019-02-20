import numpy as N
from collections import namedtuple
import os

class Furby_Error(Exception):
  def __init__(self, message, id):
    self.id = id
    self.message = message

class Furby_reader:
    def __init__(self, filename):
        self.filename=filename
        self.read_header(filename)

    def read_header(self, filename):
      if not os.path.exists(filename):
        raise Furby_Error("The furby file ({0}) does not exist".format(filename), 0)
      fil = open(filename, 'r')

      x=fil.read(16384).strip('\x00').strip()
      if "ID" not in x:
        raise Furby_Error("File: '{0}' does not seem to have a valid furby header".format(self.filename), 0)
      h=x.split("\n")
      Header_tmp={}

      for line in h:
          key=line.split()[0]
          val=line.split()[1]
          cval=self.check_type(val)
          if key == 'ID':
            cval = str(val)
          Header_tmp[key] = cval

      keys=Header_tmp.keys()
      values=Header_tmp.values()
      tmp=namedtuple("HEADER", keys)
      self.header=tmp(*values)


    def read_data(self, start=0, dd=False):
	#assume we always need to read the complete furby, i.e. nsamps = -1
        fil=open(self.filename)
        fil.seek(self.header.HDR_SIZE + N.max([0, start]))

        count=-1
        allowed_nbits = N.array([32, 16, 8])
	dtypes = N.array(['float32', 'uint16', 'uint8'])	#Assume 32 bit is always float32, and not uint32

	dtype = dtypes[  N.where(allowed_nbits == self.header.NBIT)[0][0]   ]
        data=N.fromfile(fil, count=count, dtype=dtype)
        fil.close()
	self.data = self.reshape_data(data)
	if dd:
	  self.data = self.dedisperse(self.data)
	return self.data

    def reshape_data(self, data):
	if self.header.ORDER == "TF":
	    d = data.reshape(-1, self.header.NCHAN).T
	elif self.header.ORDER == "FT":
	    d = data.reshape(self.header.NCHAN, -1)
	else:
	    raise Furby_Error("Unsupported ORDER in input : {0}".format(self.header.ORDER), 0)
	return d
	    
    def dedisperse(self, data, dm = None):
	if not dm:
	    dm = self.header.DM
	chw = self.header.BW / self.header.NCHAN
	foff = chw / 2.
	if self.header.BW>0:
	    f0 = self.header.FBOTTOM
	if self.header.BW<0:
	    f0 = self.header.FTOP

	tsamp = self.header.TSAMP / 1e6			#Dada header has to have tsamp in usec

	fch = (f0 + foff) + N.arange(self.header.NCHAN) * chw	#(f0 + foff) becomes the centre frequency of the first channel
	delays = dm * 4.14881e3 * ( fch**-2 - (f0+foff)**-2 )	#in seconds
	delays -= delays[int(self.header.NCHAN/2)]
	delays_in_samples = N.rint(delays / tsamp).astype('int')

	d_data = []
	for i, row in enumerate(data):
	    d_data.append(N.roll(row, -1*delays_in_samples[i]))
	d_data = N.array(d_data)
	return d_data	


    def check_type(self, val):
        try:
            ans=int(val)
            return ans
        except ValueError:
            try:
                ans=float(val)
                return ans
            except ValueError:
                if val.lower()=="false":
                    return False
                if val.lower()=="true":
                    return True
                else:
                    return val

