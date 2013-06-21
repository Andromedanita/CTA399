from manage import *

data=LoadData('ForMartenCleaned350Trim')
del data[:15]
mjd=IterativeFloatAppend(data,3)
int_mjd=[]
dec_mjd=[]
#for i in range(len(mjd)):
    #result=mjd[i].partition(".")
    #mjd_int=result[0]
    #int_mjd.append(int(mjd_int))
    #mjd_dec=result[2]
    #dec_mjd.append(float(mjd_dec))

t=Time(mjd,format='mjd',scale='tt',precision=9)
output=t.iso

WriteFile(output,'350Trim.dat')
