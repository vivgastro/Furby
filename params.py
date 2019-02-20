#Config file containing the parameters to be used in furby_generator

#------ Instrument params ----------

ftop                  1500.0      #MHz
fbottom               1200.0      #MHz
nch                   400
#tsamp                 0.00032768  #s
tsamp                 0.00065536  #s
noise_per_channel     1.0


#------ Astrophysical Constants ------
D                     4.14881e6   #ms ; from Pulsar Handbook Section 4.1.1
scattering_index      4.4


#------ Furby specific params --------
N_per_IDseries        1000
ID_series             1


