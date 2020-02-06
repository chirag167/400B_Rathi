#!/usr/bin/env python
# coding: utf-8

# In[2]:


# importing the numpy and astropy modules in lines 2-4.
import numpy as np

from astropy import units as u


def Read(filename): # defining a function that will read filename.
    file = open(filename,'r') # opens filename and reads it ('r')
    line1 = file.readline() # reads line 1.
    label, value = line1.split() # splits line 1 at whitespaces. 
    time = float(value)*u.Myr # stores value in a variable time.
    
    line2 = file.readline() # similarly, reads line 2
    label, value = line2.split() # splits line 2 at whitespaces.
    particles = float(value) # stores value in a variable particles.
    
    file.close() # close the file.
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header = 3) # save the remainder of the file as arrays (under the variable, data).
    
    return time,particles,data # return the variables, time, particles and data.
    


time,particles,data = Read("MW_000.txt") # calling the Read function to read "MW_000.txt" and store the ouput in time,particles and data.


# In[ ]:




