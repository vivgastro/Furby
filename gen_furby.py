#!/home/observer/miniconda2/bin/python

import numpy as N
import matplotlib.pyplot as M
import argparse
import os, sys
import glob
import time

from collections import namedtuple
from parse_cfg import parse_cfg as pcfg

P = pcfg("params.cfg")

consts = {
    'tfactor': int( P.tsamp / 0.00001024 ),               #40.96 microseconds
    'ffactor': int( ((P.ftop-P.fbottom)/P.nch)/0.01),      #We dont want dmsmear to be approximated beyond 5 kHz chw. So ffactor = chw/ 0.005
    }

tmp = namedtuple("co", consts.keys())
C = tmp(*consts.values())


def tscrunch(data, tx):
    if tx==1:
        return data
    
    if len(data.shape)==1:
        endpoint = int(len(data) / tx) * tx
        return data[:endpoint].reshape(-1, tx).sum(axis=-1)

    if len(data.shape)==2:
        nr=data.shape[0]
        nc=data.shape[1]
    
        endpoint=int(nc/tx) * tx
        tmp=data[:,:endpoint].reshape(nr, nc/tx, tx)
        tsdata=tmp.sum(axis=-1)

        return tsdata
    else:
        raise RuntimeError("Can only scrunch 1D/2D arrays")

def fscrunch(data, fx):
    if fx==1:
      return data
    if fx==0:
      raise ValueError("Cannot fscrunch by a factor of 0")
    nr = data.shape[0]
    nc = data.shape[1]

    if nr%fx!=0:
      raise RuntimeError("Cannot scrunch at factors which do not exactly divide the no. of channels")
    fsdata = N.mean(data.reshape(nr/fx, -1, nc), axis=1)
    return fsdata

def get_clean_noise_rms():
    #noise rms per channel
    return (P.noise_per_channel)

def gauss(x, a, x0, sigma):
    return a/N.sqrt(2*N.pi*sigma**2) * N.exp(-(x-x0*1.)**2 / (2.*sigma**2))

def gauss2(x, a, x0, FWHM):
    sigma = FWHM/2. /(2*N.log(2))**0.5			#FWHM = 2 * sqrt( 2 * ln(2) ) * sigma
    return a/N.sqrt(2*N.pi*sigma**2) * N.exp(-(x-x0*1.)**2 / (2.*sigma**2))


def get_pure_frb(snr, width, nch, nsamps):
    clean_noise_rms = get_clean_noise_rms()     #this is ideal noise rms per channel
    snr_per_channel = snr*1./N.sqrt(nch)        #Dividing snr equally among all channels for the pure case

    tmp_sigma = width/2. /(2*N.log(2))**0.5	#width is supposed to be FWHM
    W_tophat_gauss = N.sqrt(2*N.pi) * tmp_sigma

    desired_signal = snr_per_channel * clean_noise_rms * N.sqrt(W_tophat_gauss)

    x=N.arange(int(nsamps * C.tfactor) )
    width = width * C.tfactor
    pure_frb_single_channel = gauss2(x, desired_signal, int(len(x)/2), width)

    if N.abs(N.sum(pure_frb_single_channel) - desired_signal) > desired_signal/50.:
      raise RuntimeError("The generated signal is off by more than 2% of the desired value, desired_signal = {0}, generated_signal = {1}. Diff: {2}%".format(desired_signal, N.sum(pure_frb_single_channel), ((N.sum(pure_frb_single_channel) - desired_signal)/desired_signal * 100) ))

    pure_frb = N.array([pure_frb_single_channel] * nch)     #Copying single channel nch times as a 2D array
    
    assert pure_frb.shape[0] == nch, "Could not copy 1D array {0} times".format(nch)
    
    return pure_frb

def get_bandpass(nch):
    bp = N.loadtxt("/home/vgupta/resources/BANDPASS_normalized_320chan.cfg")
    if nch == 320:
        pass
    elif nch == 40:
        bp = tscrunch(bp, 8) / 8.
    else:
        raise ValueError("NCHAN expected: [40 or 320]. Got: {0}".format(str(nch)))
    return bp*1./bp.max()

def apply_bandpass(frb):
    nch = frb.shape[0]
    bp = get_bandpass(nch)
    bp = bp/bp.max()
    #bp = bp - bp.mean() +1

    frb = frb * bp.reshape(-1,1)
    return frb

def create_freq_structure(frb, kind):
    nch = frb.shape[0]
    x = N.arange(nch)
    #kind of scintillation
    
    if kind == 'flat':
        f = N.ones(nch)
    if kind == 'slope':
        slope = N.random.uniform(-0.5*nch, 0.5*nch, 1)  #Slope will be a random number between -0.5 and 0.5
        f = x * slope
    if kind == 'smooth_envelope':
        center = N.random.uniform(0, nch, 1)            #Location of Maxima of the smooth envelope can be on any channel
        z1 = center - nch/2
        z2 = center + nch/2
        f = -1 * (x - z1) * (x - z2)
    if kind == 'two_peaks':
        z1 = 0
        z2 = N.random.uniform(0 + 1, nch/2, 1)
        z3 = N.random.uniform(nch/2, nch-1 , 1)
        z4 = nch
        f = -1 * (x-z1) * (x-z2) * (x-z3) * (x-z4)
    if kind == 'three_peaks':
        z1 = 0
        z2 = N.random.uniform(0 +1, nch/4, 1)
        z3 = N.random.uniform(1*nch/4, 2*nch/4, 1)
        z4 = N.random.uniform(2*nch/4, 3*nch/4, 1)
        z5 = N.random.uniform(3*nch/4, nch-1, 1)
        z6 = nch
        f = -1 * (x-z1) * (x-z2) * (x-z3) * (x-z4) * (x-z5) * (x-z6)
    if kind == 'ASKAP':
        n_blobs = N.floor(N.random.exponential(scale = 3, size=1)) + 1
        f = N.zeros(nch)
        for i in range(n_blobs):
            center_of_blob = N.random.uniform(0, nch, 1)
            #We want roughly 4 +- 1 MHz blobs. 4 MHz = 4/chw chans = 4./((P.ftop - P.bottom)/nch) chans
            NCHAN_PER_MHz = N.abs(1./( (P.ftop-P.fbottom)/nch ))
            width_of_blob = N.random.normal(loc = 4.*NCHAN_PER_MHz, scale = NCHAN_PER_MHz, size = 1)  
            amp_of_blob = N.random.uniform(1, 3, 1)             #For just one blob (n_blobs=1), this does not matter because we rescale the maxima to 1 evetually. For more than one blobs, this random amp will set the relative power in different blobs. So, the power in weakest blob can be as low as 1/3rd of the strongest blob)
            f += gauss(x, amp_of_blob, center_of_blob, width_of_blob)
            
    if kind != 'flat':
      f = f - f.min()                                     #Bringing the minima to 0
      f = f * 1./f.max()                                  #Bringing the maxima to 1
      f = f - f.mean() + 1                                #Shifting the mean to 1

    frb = frb * f.reshape(-1, 1)
    return frb, f

def disperse(frb, dm, pre_shift, dmsmear):
    tsum = 0
    if not dmsmear:
      ffactor = 1
    else:
      ffactor = C.ffactor

    if args.v:
      print "Ffactor:", ffactor
    nch = frb.shape[0] * ffactor
    tres = P.tsamp / C.tfactor *1e3  #ms. Effectively 10.24 micro-seconds, just framing it in terms of tres of Hires filterbanks

    chw = (P.ftop-P.fbottom)/nch
    f_ch = N.linspace(P.ftop - chw/2, P.fbottom + chw/2, nch)
    delays = P.D * f_ch**(-2) * dm    #Only works if freq in MHz and D in ms. Output delays in ms
    delays -= delays[int(nch/2)]
    delays_in_samples = N.rint(delays / tres).astype('int') #here we will have slight approximations due to quantization, but at 10.24 usec resolution they should be minimal

    #nsamps = delays_in_samples[-1] - delays_in_samples[0] + 2*frb.shape[1]
    #nsamps = delays_in_samples[-1]*2 + 2*frb.shape[1]
    nsamps = 9000 * C.tfactor
    start = nsamps/2 - int(pre_shift*C.tfactor)
    end = start + frb.shape[1]

    dispersed_frb = N.zeros(nch * nsamps).reshape(nch, nsamps)
    undispersed_time_series = N.zeros(nsamps)
    idxs = N.arange(nsamps)

    if args.v:
      print "Initial frb shape", frb.shape, "nsamps:",nsamps, "\nstart, end",start, end, "Tfactor, Ffactor, pre_shift", C.tfactor, ffactor, int(pre_shift)

    for i in range(nch):
      delay = delays_in_samples[i]
      mid_channel_coarse = int(i/ffactor) *ffactor + int(ffactor/2.0)
      dispersed_frb[i, start+delay: end+delay] += frb[int(i/ffactor)]
      undispersed_time_series += N.take(dispersed_frb[i], idxs + delays_in_samples[mid_channel_coarse], mode='wrap')
      tsum = undispersed_time_series.sum()
      if args.v:
        sys.stdout.write("nch : {0}/{1}  tsum = {2}\r".format(i, nch, tsum))
    final_dispersed_frb = fscrunch(dispersed_frb, ffactor)
    return final_dispersed_frb, undispersed_time_series/ffactor, nsamps/C.tfactor


def scatter(frb, tau0, nsamps, desired_snr): 
    nch = frb.shape[0]
    ftop = P.ftop        #MHz
    fbottom = P.fbottom     #MHz
    chw = (ftop-fbottom)/(1.0*nch)
    f_ch = N.linspace(ftop - chw/2, fbottom + chw/2, nch)
    nsamps = nsamps * C.tfactor
    tau0 = tau0 * C.tfactor

    k = tau0 * (f_ch[0])**P.scattering_index     #proportionality constant
    taus = k / f_ch**P.scattering_index          #Calculating tau for each channel
    exps=[]
    scattered_frb=[]
    for i,t in enumerate(taus):
        exps.append( N.exp(-1 * N.arange(nsamps) / t) )       #making the exponential with which to convolve each channel
        result = N.convolve(frb[i], exps[-1])                 #convolving each channel with the corresponding exponential ( N.convolve gives the output with length = len(frb) + len(exp) )
        #result *= 1./result.max() * frb[i].max()
        result *= 1./result.sum() * frb[i].sum()
        scattered_frb.append(result) 

    scattered_frb=N.array(scattered_frb)
    scattered_tseries = scattered_frb.sum(axis=0)
    scattered_width = scattered_tseries.sum() / N.max(scattered_tseries) / C.tfactor
    new_snr = scattered_tseries.sum() / (N.sqrt(nch) * get_clean_noise_rms() ) / N.sqrt(scattered_width)
    normalizing_factor = new_snr / desired_snr
    scattered_frb /= normalizing_factor
    return scattered_frb

def make_psrdada_header(hdr_len, params):
    header=""
    for i in params:
        header += i
        tabs = 3 - int(len(i)/8)
        header += "\t"*tabs
        header += str(params[i])+"\n"
    leftover_space = hdr_len - len(header)
    header += '\0' * leftover_space
    return header

def get_FWHM(frb_tseries):
    maxx = N.argmax(frb_tseries)
    hp = frb_tseries[maxx] / 2.
    #Finding the half-power points
    hpp1 = (N.abs(frb_tseries[:maxx] - hp)).argmin()
    hpp2 = (N.abs(frb_tseries[maxx:] - hp)).argmin() + maxx

    FWHM = (1.*hpp2-hpp1)/C.tfactor
    assert FWHM>0, "FWHM calculation went wrong somewhere. HPP points, maxx point and FWHM are {0} {1} {2} {3}".format(hpp1, hpp2, maxx, FWHM)
    return FWHM, N.max(tscrunch(frb_tseries, C.tfactor))

def start_logging(ctl, db_d):
    if os.path.exists(ctl):
	logger = open(ctl, 'a')
    else:
	logger = open(ctl, 'a')
	logger.write("#This is the furby catalogue for {0} directory\n".format(db_d))
	logger.write("#Created on : {0}\n\n".format(time.asctime()))
	existing_furbies = glob.glob(db_d+"furby_*")
	if len(existing_furbies) > 0:
	    logger.write("#The following furbies were found to be present in the directory before the creation of this catalogue:\n")
	    for furby in existing_furbies:
		logger.write("#"+furby+"\n")
	    logger.write("\n")
	logger.write("#FURBY_ID\tDM\tFWHM\tTAU0\tSNR\tSIGNAL\n")
    return logger

def check_for_permissions(db_d):
    if not os.path.exists(db_d):
	try:
	    print "The database directory: {0} does not exist. Attempting to create it now.".format(db_d)
	    os.makedirs(db_d)
	except OSError as E:
	    print "Attempt to create the database directory failed because:\n{0}".format(E.strerror)
	    print "Exiting...."
	    sys.exit(1)

    if os.access(db_d, os.W_OK) and os.access(db_d, os.X_OK):
        return 
    else:
	print "Do not have permissions to write/create in the database directory: {0}".format(db_d)
	print "Exiting..."
	sys.exit(1)

def main(args):
    database_directory = args.D
    if not database_directory.endswith("/"):
	database_directory += "/"
    
    order = "TF"
    
    if args.plot and args.Num > 1:
      raise IOError("Sorry cannot plot more than one furby at a time")

    if not args.plot:
	check_for_permissions(database_directory)
	catalogue = database_directory+"furbies.cat"
	logger = start_logging(catalogue, database_directory)

    ID_series = P.ID_series

    if args.v:
        print "Starting FRB Generator..."
    tsamp = P.tsamp                              #seconds
    nch = P.nch
    supported_kinds = ["flat", "slope", "smooth_envelope", "two_peaks", "three_peaks", "ASKAP"]
    
    if isinstance(args.snr, float):
      snrs = args.snr * N.ones(args.Num)
    elif isinstance(args.snr, list) and len(args.snr) ==1:
      snrs = args.snr[0] * N.ones(args.Num)
    elif isinstance(args.snr, list) and len(args.snr) ==2:
      snrs = N.random.uniform(args.snr[0], args.snr[1], args.Num)
    else:
      raise IOError("Invalid input for SNR")
    #snr = 15

    if isinstance(args.width, float):
      #widths = (args.width *1e-3/ P.tsamp) * N.ones(args.Num)
      widths = args.width *1e-3/ P.tsamp * N.ones(args.Num)
    elif isinstance(args.width, list) and len(args.width) ==1:
      #widths = int(args.width[0]*1e-3/P.tsamp) * N.ones(args.Num)
      widths = args.width[0]*1e-3/P.tsamp * N.ones(args.Num)
    elif isinstance(args.width, list) and len(args.width) ==2:
      #widths = N.random.randint(args.width[0]*1e-3/P.tsamp,args.width[1]*1e-3/P.tsamp, args.Num)   #samples
      widths = N.random.uniform(args.width[0]*1e-3/P.tsamp,args.width[1]*1e-3/P.tsamp, args.Num)   #samples
    else:
      raise IOError("Invalid input for Width")
    #width =3
    
    if isinstance(args.dm, float):
      dms = args.dm * N.ones(args.Num)
    elif isinstance(args.dm, list) and len(args.dm) ==1:
      dms = args.dm[0] * N.ones(args.Num)
    elif isinstance(args.dm, list) and len(args.dm) ==2:
      dms = N.random.uniform(args.dm[0], args.dm[1], args.Num)
    else:
      raise IOError("Invalid input for DM")
    #dm = 900
# -----------------------------------------------------------------------------------------------------------------
    for num in range(args.Num):
      dm = dms[num]
      snr = snrs[num]
      width = widths[num]

      if args.kind:
          kind = args.kind
          assert kind in supported_kinds
      else:
        kind = supported_kinds[N.random.randint(0, len(supported_kinds), 1)[0]]

      while(True):
          ID = N.random.randint(( ID_series*P.N_per_IDseries + 1), (ID_series+1)*P.N_per_IDseries, 1)[0] 
          ID = str(ID).zfill(int(N.log10(P.N_per_IDseries)))
          name = "furby_"+ID
          if os.path.exists(database_directory+name):
              continue
          else:
              break
      if args.scatter:
        tau0 = N.abs(N.random.normal(loc = dm / 1000., scale = 2, size=1))[0] 
      else:
        tau0 = 0
      #tau0 = 10.1/C.tfactor
      
      nsamps_for_gaussian = 5 * width            # = half of nsamps required for the gaussian. i.e. The total nsamps will be 10 * sigma.
      if nsamps_for_gaussian < 1:
        nsamps_for_gaussian = 1
      nsamps_for_exponential = int(6 * tau0 * ((P.ftop+P.fbottom)/2 / P.fbottom)**P.scattering_index)

      if args.v:
          print "Randomly generated Parameters:"
          print "ID= {0}, SNR= {1}, Width= {2}ms, DM= {3}, tau0= {4}ms, kind= {5}".format(ID, snr,width*tsamp*1e3, dm, tau0*tsamp*1e3, kind)
      
      if args.v:
          print "Getting pure FRB"
      try:
        frb = get_pure_frb(snr=snr, width = width, nch=nch, nsamps=nsamps_for_gaussian)
      except RuntimeError as R:
        print R
        continue
      
      pure_signal = N.sum(frb.flatten())
      if args.v:
          print "Creating frequency structure"
      frb,f = create_freq_structure(frb, kind=kind)

      pure_width = pure_signal / N.max( frb.sum(axis=0) )   /C.tfactor
      pure_snr = pure_signal / (N.sqrt(nch) * get_clean_noise_rms() * N.sqrt(pure_width))
      
      if args.v:
        print "Pure signal (input) = {0}, signal after freq_struct = {1}, pure_snr = {2}, pure_width = {3}ms".format(pure_signal, N.sum(frb.flatten()), pure_snr, pure_width * P.tsamp * 1e3)  

      #if args.v:
      #  print "Applying Bandpass"
      #frb = apply_bandpass(frb)
      #if args.v:
      #  print "Signal after bandpass calib = {0}".format(N.sum(frb.flatten()))
      
      if args.v:
          print "Scattering the FRB"
      if nsamps_for_exponential==0:
        print "Tau0 = 0, no scattering applied"
        pass
      else:
        frb = scatter(frb, tau0, nsamps_for_exponential, pure_snr)
      sky_signal = N.sum(frb.flatten())
      sky_frb_tseries = N.sum(frb, axis=0)
      #sky_frb_peak = N.max( tscrunch(sky_frb_tseries, C.tfactor)   )
      sky_frb_top_hat_width = sky_signal / N.max(sky_frb_tseries) / C.tfactor
      #sky_frb_top_hat_width = sky_signal / sky_frb_peak
      sky_snr = sky_signal / ( get_clean_noise_rms() * N.sqrt(nch) * N.sqrt(sky_frb_top_hat_width) )

      if args.v:
        print "Sky_signal = {0}, sky_width = {1}ms, sky_snr = {2}".format(sky_signal, sky_frb_top_hat_width * P.tsamp * 1e3, sky_snr)
      
      frb_b_d = frb.copy()      #FRB before dispersing
      
      if args.v:
          print "Dispersing"
      frb, undispersed_tseries, NSAMPS = disperse(frb, dm, pre_shift = nsamps_for_gaussian,dmsmear = args.dmsmear) #Remember, nsamps_for_gaussian is already half the number of samples
      signal_after_disp = N.sum(frb.flatten())

      FWHM,maxima = get_FWHM(undispersed_tseries)

      if args.v:
          print "Scrunching"
      scrunched_frb = tscrunch(frb, C.tfactor)
      signal_after_scrunching = N.sum(scrunched_frb.flatten())

      #final_top_hat_width = signal_after_scrunching / maxima		#comes out in samples
      final_top_hat_width = signal_after_scrunching / N.max(undispersed_tseries) / C.tfactor
      output_snr = signal_after_scrunching / (get_clean_noise_rms() * N.sqrt(nch) *  N.sqrt(final_top_hat_width) )

      if args.v:
          print "Input signal, Sky_signal, Output signal, Input SNR, Sky SNR, Output SNR, Final_top_hat_width",  pure_signal, sky_signal, signal_after_scrunching, snr, sky_snr, output_snr, final_top_hat_width * tsamp * 1e3, "ms\n"
      
      final_frb = scrunched_frb.astype('float32')

      if args.v:
          print "Generating Header"
      header_size = 16384                 #bytes
      params = {
              "HDR_VERSION":  1.0,
              "HDR_SIZE": header_size,
              "TELESCOPE":    "MOST",
              "ID":   ID,
              "SOURCE":	name,
              "FREQ": (P.ftop + P.fbottom)/2,         #MHz
              "BW":   P.fbottom - P.ftop,              #MHz. Negative to indicate that the first channel has highest frequency
              "NPOL": 1,
              "NBIT": 32,
              "TSAMP":    tsamp * 1e6,      	#usec	--Dont change this, dspsr needs tsamp in microseconds to work properly
              "NSAMPS":	NSAMPS,
              "UTC_START":    "2018-02-12-00:00:00",       #Never happened :P
              "STATE":    "Intensity",
              "OBS_OFFSET":   0,           #samples
              "NCHAN":    nch,
              "ORDER":    order,            #This is ensured later in the code while writing data with flatten(order=O)
              "FTOP": P.ftop,        #MHz
              "FBOTTOM": P.fbottom,     #MHz
              "TRACKING": "false",
              "FB_ENABLED":   "true",
              "INSTRUMENT":   "MOPSR",
              "SNR":  output_snr,
              "SKY_SNR":  sky_snr,
      	      "WIDTH":	final_top_hat_width*tsamp*1e3,	#milliseconds
              "SIGNAL":   signal_after_scrunching,
              "SKY_SIGNAL":   sky_signal,
              "WIDTH_GAUSS":  width*tsamp*1e3,             #milliseconds
              "FWHM": FWHM*tsamp*1e3,                      #milliseconds
              "DM":   dm,                              #pc / cm^3
              "DMSMEAR":  str(args.dmsmear).upper(),
              "TAU0": tau0*tsamp*1e3,                      #milliseconds
              "KIND": kind,                            #kind of frequency structure put into the FRB
              "PEAK": maxima,                          #peak signal of the FRB
              }
      header_string = make_psrdada_header(header_size, params)
      if not args.plot:
          if args.v:
              print "Saving the FRB in '"+database_directory+"' directory"

          out = open(database_directory+name, 'wb')
          out.write(header_string)
          if order == "TF":
              O = 'F'				#Order = 'F' means column-major, since we are writing data in TF format.
          if order == "FT":
              O = 'C'				#Order = 'C' means row-major, since we are writing data in FT format.
          final_frb.flatten(order=O).tofile(out)          
          out.close()
          logger.write(ID+"\t"+str(dm)+"\t"+str(FWHM*tsamp)+"\t" + str(tau0*tsamp)+"\t" + str(output_snr)+"\t"+str(signal_after_scrunching)+"\n")
          
          print "Name : ", params["SOURCE"]
#------------------------------------------------------------------------------------------------------------
    if not args.plot:
        logger.close()

    if args.plot:
        if args.v:
            print "Plotting"
        M.figure()
        M.imshow(frb_b_d, aspect='auto', cmap='afmhot', interpolation='None')
        M.title("FRB before dispersion")

        M.figure()
        M.plot(N.sum(frb_b_d, axis=0))
        M.title("FRB tseries before dispersing")

        M.figure()
        M.imshow(scrunched_frb, aspect = 'auto', cmap = 'afmhot', interpolation='nearest')
        M.title("FRB")

        M.figure()
        M.plot(tscrunch(undispersed_tseries, C.tfactor))
        M.title("Un-dispersed FRB time series. We should be aiming to recover this profile")
        M.show()

if __name__ == '__main__':
    a=argparse.ArgumentParser()
    a.add_argument("Num", type=int, help="Number of furbies to generate")
    a.add_argument("-kind", type=str, help="Kind of frequency structure wanted. Options:[slope, smooth_envelope, two_peaks, three_peaks, ASKAP]")
    a.add_argument("-plot", action='store_true', help = "Plot the FRB instead of saving it?", default = False)
    a.add_argument("-dm", nargs='+', type=float, help="DM or DM range endpoints", default = 1000.0)
    a.add_argument("-snr", nargs='+', type=float, help="SNR or SNR range endpoints", default = 20.0)
    a.add_argument("-width", nargs='+', type=float, help="Width or width range endpoints (in ms)", default = 2.0)
    a.add_argument("-dmsmear", action='store_true', help = "Enable smearing within individual channels (def=False)", default=False)
    a.add_argument("-scatter", action='store_true', help='Enable scattering (def = False)', default=False)
    a.add_argument("-v", action='store_true', help="Verbose output", default = False)
    a.add_argument("-D", type=str, help="Path to the database to which the furby should be added (def = ./)", default = "./")
    args= a.parse_args()
    main(args)
