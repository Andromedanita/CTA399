import sys
import numpy as np
import matplotlib as plt
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.mlab as mlab
import scipy
from scipy.optimize import leastsq
import pylab

#This program gets the values of apparent major axis, inclination, time, period and longitude of ascending node as input and generates the elliptical orbit of the pulsar. We assumed that the true orbit is circular but it would be elliptical the way we see it unless inclination is 90 degrees. After that, it will randomize the values of right ascention and declination and finally, generates a fit to the data.

min_x=-90
max_x=90
step=0.8
x=mgrid [min_x: (max_x+step): step]

#Inputs user must enter
a = int(sys.argv[1]) #projected major axis
i = int(sys.argv[2]) #inclination
t= int(sys.argv[3])  #time 
T= int(sys.argv[4])  #period of the pulsar orbit
w= int(sys.argv[5])  #longitude of asceding node

t=linspace(-3*T,3*T,20)

#minor axis of the ellipse
b=a*(np.cos(i))

#equations for inclination
x1=a*(np.cos(t/T))
y1=b*(np.sin(t/T))

#equations for longitude of ascending node(rotation matrix transformation)
x=(x1*np.cos(w))-(y1*np.sin(w))
y=(x1*np.sin(w))-(y1*np.cos(w))

#plotting the orbit
plt.plot(x,y, '*', color='r')

#Using scipy.optimize.leastsquare for question 3

#defining error for x values
def residual1(p,x,t):
    i,w = p
    b=a*(np.cos(i))
    x1=a*(np.cos(t/T))
    y1=b*(np.sin(t/T))
    err=x-(x1*np.cos(w))-(y1*np.sin(w))
    return err

#defining error for y values
def residual2(p,y,t):
    i,w = p
    b=a*(np.cos(i))
    x1=a*(np.cos(t/T))
    y1=b*(np.sin(t/T))
    err=y-(x1*np.sin(w))-(y1*np.cos(w))
    return err

#p0 is the initial guess
p0=[i,w]
out_x=leastsq(residual1, p0, args=(x, t))
out_y=leastsq(residual2, p0, args=(y, t))
print out_x[0][0]

#The fit ellipse equations
#b_f is the minor axis of the fit ellipse
b_f=a*(np.cos(out_x[0][0]))

#equation of fit ellipse only considering inclination
x2=a*(np.cos(t/T))
y2=b_f*(np.sin(t/T))

#equation of the fit ellipse with error values and considering longitude of 
#ascending node
x_f=(x2*np.cos(out_y[0][1]))-(y2*np.sin(out_y[0][1]))
y_f=(x2*np.sin(out_y[0][1]))-(y2*np.cos(out_y[0][1]))
print (out_x), (out_y)

#plotting the fit ellipse
plt.plot(x_f,y_f, color='b')


#plot customization   
plt.xlabel("Declination(deg)")
plt.ylabel("Right Ascention(deg)")
plt.title("Position of the pulsar")
#chose these values ti show declination and right ascension
plt.axis([-90, 90, -180, 180])
plt.legend(('Actual ellipse','fit ellipse') , loc='upper right')
plt.locator_params(nbins=10)
#plt.text(i, w, fontsize=20, fontname='Times New Roman')
print "i=", (i), "w=", (w)
plt.show()
