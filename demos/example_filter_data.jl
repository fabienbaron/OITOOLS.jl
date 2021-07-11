#
# OIFITS filter example
#
using OITOOLS

# read the data file
oifitsfile = "./data/AlphaCenA.oifits";


data = (readoifits(oifitsfile))[1,1]
uvplot(data, bywavelength=true)
# Create the list of uv points and corresponding observables to remove
filter = set_data_filter(data, filter_wav=true, minwav=1.5e-6, maxwav=1.68e-6)
data2 = filter_data(data, filter) # data2 is filtered
uvplot(data2, bywavelength=true)

# One can use multiple filters at the same time
filter = set_data_filter(data, filter_wav=true, minwav=1.5e-6, maxwav=1.68e-6, filter_mjd = true, minmjd=57532.07, maxmjd = 57532.3 )
data3 = filter_data(data, filter)

# "Bad" data (NaN, negative errors, etc.) and flagged data can be read in first
# then filtered out afterwards
data = (readoifits(oifitsfile,filter_bad_data=false))[1,1]
filter = set_data_filter(data, filter_bad_data=true)
data3 = filter_data(data, filter)
