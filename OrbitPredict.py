from scipy.optimize import leastsq
import numpy as np
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt
import argparse

#Create an command line argument parser. 
parser=argparse.ArgumentParser()

#Add command line options
#Required arguments
parser.add_argument('--T', type=float, required=True, help=''' The period of the orbit in seconds. ''')
parser.add_argument('--a', type=float, required=True, help=''' The length of the projected semi-major axis of the system in arcseconds.''')
parser.add_argument('--data',type=file, required=True, help=''' A file containing data points to be fitted. First column should be right ascenion, next column the errors therein, the third declination, and the fourth errors therein.''')

#Optional arguments
#time of ascending node defaults to zero for the sake of simplicity
parser.add_argument('--t_0', type=float, default=0, help=''' The time of ascending node as a Julian Date.''') 
#proper motion in either direction - defaults to zero
parser.add_argument('--rapm', type=float, default=0, help=''' The proper motion of the orbit in the direction of right ascension in arcseconds per second.''')
parser.add_argument('--decpm', type=float, default=0, help=''' The proper motion of the orbit in the direction of declination in arcseconds per second.''')


#Parse the command line arguments
args=parser.parse_args()

#Defining variables from input
P=args.T #period of orbit
maj=args.a #projected semi major axis
t_0=args.t_0 #time of ascending node
ra_pm=args.rapm #proper motion in the direction of right ascension
dec_pm=args.decpm #proper motion in the direction of declination
data=args.data #represent input file as string

#Need right ascension, declination and their errors from input file

contents=np.loadtxt(data)
ra=contents[:,0]
ra_err=contents[:,1]
dec=contents[:,2]
dec_err=contents[:,3]

time=np.linspace(t_0,(orb*2*P*np.pi)+t_0,len(ra))

plt.errorbar(ra,dec,xerr=ra_err,yerr=dec_err,fmt='bo',linestyle='none', label='Measured Data')


#DEFINING FUNCTIONS

#Define a function that take time, initial time, semi-major axis, period, proper
#motion in the directions of right ascension and declination, inclination and 
#longitude of ascending node and creates a list of points in the orbit as it 
#appears on the sky

def apparent_orbit(t,ti,a,p,v,u,ang,rot):
    #Define the minor axis of the ellipse in terms of the major axis and the 
    #inclination 
    b=a*np.cos(ang)
    #Write basic parametric equations for an ellipse
    x=a*np.cos((t-ti)/p) 
    y=b*np.sin((t-ti)/p)
    #Now rotate the ellipse according to the longitude of the ascending node 
    #and incorporate the proper motion
    mod_x=-x*np.sin(rot)+y*np.cos(rot)+v*(t-ti)
    mod_y=x*np.cos(rot)+y*np.sin(rot)+u*(t-ti)
    return (mod_x,mod_y)

#Choose random modified x and y from a Gaussian centred at the original 
#point with an standard deviation specified by 'scale' - essentially, choose
#a point within the errors of the original one.

def Randomize(x,y,x_err,y_err):
    mod_x=np.random.normal(loc=x,scale=x_err)
    mod_y=np.random.normal(loc=y,scale=y_err)
    return (mod_x,mod_y)


#A function that calculates the difference between input y values and the y 
#values of an ellipse with parameters ang and rot, two angles that describe its
#eccentricity and rotation around the origin, respectively
#residuals_y takes two parameters, y-values, and time and returns a list
    
def residuals_y(p,y,t):
    #ang and rot are two angles described in the parameter list p
    ang,rot=p
    #the minor axis of the ellipse is related to the major axis by the following
    mi=maj*np.cos(ang)
    #basic parametric equations for an ellipse
    el_x=maj*np.cos((t-t_0)/P)
    el_y=mi*np.sin((t-t_0)/P)
    #modfied y value that has been rotated around the origin and had proper 
    #motion incorporated
    proper_dist=(dec_pm*(t-t_0))
    mod_y=el_x*np.cos(rot)+el_y*np.sin(rot)+proper_dist
    #the difference between the input y values and y values on an ellipse
    err=y-mod_y
    return err

#A function that calculates the difference between input x values and the x 
#values of an ellipse with parameters ang and rot, two angles that describe its
#eccentricity and rotation around the origin, respectively
#residuals_x takes two parameters, x-values and time and returns a list
    
def residuals_x(p,x,t):
    #ang and rot are two angles described in the parameter list p
    ang,rot=p
    #the minor axis of the ellipse is related to the major axis by the following
    mi=maj*np.cos(ang)
    #basic parametric equations for an ellipse
    el_x=maj*np.cos((t-t_0)/P)
    el_y=mi*np.sin((t-t_0)/P)
    #modified x value that has been rotated around the origin and had proper 
    #motion incorporated
    proper_dist=(ra_pm*(t-t_0))
    mod_x=-el_x*np.sin(rot)+el_y*np.cos(rot)+proper_dist
    #the difference between the input x values and the x values on an ellipse
    err=x-mod_x
    return err


#RANDOMIZE THE VALUES IN THE ELLIPSE

j=0 #an index representing the number of times the fitter has run
parameters=[] #an empty list to hold the results of each iteration

#initial guess for parameters. p0[0] is inclination, p0[1] is longitude
#this guess is used by the leastsq function
#p0=[float(raw_input('''The least-squares fitting needs an intial guess to run
#0.707 is usually a good start. What is a likely inclination angle in radians?''#')),float(raw_input('''And the longitude of ascending node in radians?'''))]

p0=[0.707,0.707] #******************************************************
error_ra=np.array(ra_err)
error_dec=np.array(dec_err)

#Run the fitter 100 times
while j<100:
#Randomizing right ascension and declination

#recover a list of points from the function
    values=np.array(Randomize(ra,dec,ra_err,dec_err)) 

#slice the list into right ascension and declination
    ran_ra=values[0]#randomized right ascension
    ran_dec=values[1]#randomized declination

#FIT AN ELLIPSE TO RANDOM DATA

#minimize the squares of the y residuals, taking input y to be the randomized 
#declination and the guess as given above
#returns an array of the two predicted parameters
    dec_lsq=leastsq(residuals_y,p0,args=(ran_dec,time))

#minimize the squares of the x residuals, taking input x to be the randomized 
#right ascension and the guess as given above
#returns an array of the two predicted parameters
    ra_lsq=leastsq(residuals_x,p0,args=(ran_ra,time))

#take the average of the parameters generates by the optimizing leastsq function
    ra_err=np.sum(error_ra)/len(error_ra)
    dec_err=np.sum(error_dec)/len(error_dec)
    tot_err=ra_err+dec_err
    #want the direction with smaller error to have more influence on the 
    #estimated value. 
    ra_weight=dec_err/tot_err 
    dec_weight=ra_err/tot_err
    lsq=np.array(ra_lsq[0]*ra_weight+dec_lsq[0]*dec_weight)

#divide lsq into inclination and longitude    
    inclination=lsq[0]
    longitude=lsq[1]

#add each pair to parameters list
    values=[inclination, longitude]
    parameters.append(values)

#increment j to continue to the next iteration
    j+=1

#Extract inclination and longitude into separate lists total_inc and total_long
total_inc=[]
total_long=[]
for i in range(len(parameters)):
    total_inc.append(parameters[i][0])
    total_long.append(parameters[i][1])

#Average each list to find the values for inclination and longitude found by 
#the fitter
inc = sum(total_inc)/len(total_inc)
lon = sum(total_long)/len(total_long)

#Calculate the standard deviation over the dataset
inc_error=[]
long_error=[]
for i in range(len(parameters)):
    inc_error.append(abs(total_inc[i]-inclination))
    long_error.append(abs(total_long[i]-longitude))

inclination_error = sum(inc_error)/len(inc_error)
longitude_error = sum(long_error)/len(long_error)

#Print the results to the shell
print ''
print 'The inclination is: '
print str(inc) + '+/-' + str(inclination_error)
print ''
print 'The longitude of ascending node is: '
print str(lon) + '+/-' + str(longitude_error)

#Plot the new ellipse
time=np.linspace(t_0,(2*P*np.pi)+t_0,1000)
points=np.array(apparent_orbit(time,t_0,maj,P,ra_pm,dec_pm,inc,lon))
ra=points[0] 
dec=points[1]
plt.plot(ra,dec,color='red',label='Least-Squares Fitted Curve')
plt.xlabel('Right ascension relative to the centre of the orbit in arcseconds')
plt.ylabel('Declination relative to the centre of the orbit in arcseconds')
plt.title('Apparent orbit of a pulsar with a {0}" semimajor axis.'.format(maj))
plt.legend(('Measured Data','Least-Squares Fitted Curve'),loc='best')
plt.show()


