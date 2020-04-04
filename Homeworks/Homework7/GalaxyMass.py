#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Creating a function to calculate the total mass of a galaxy component.
# Galaxy components: Halo (DM), disk, bulge.

from ReadFile import Read
import numpy as np
from astropy import units as u

def ComponentMass(filename,particle_type):
    # Inputs:
    #      filename = name of the file to read.
    #       particle_type = galaxy component: 1 = Halo, 2 = Disk, 3 = Bulge.
    # Returns:
    #       total mass of the galaxy component.
    
    time, particles, data = Read(filename)
    # finding a particular particle_type.
    index = np.where(data['type']== particle_type)
    mass = data['m'][index]
    comp_mass = sum(mass)/1e2*u.Msun # in units of Msun/1e12. 
    
    # NOTE:In the file, the mass is already in units of Msun/1e10.
    
    return np.around(comp_mass*1e12,3) #rounding off the component mass to 3 decimal places.
    


# In[3]:


# The mass of the component 1 (Halo) for the Milky Way Galaxy is:
#HaloMass_MW = ComponentMass("MW_000.txt",1)
#print(HaloMass_MW)


# In[4]:


# The mass of component 2 (Disk stars) for the Milky Way Galaxy is:
#DiskMass_MW = ComponentMass("MW_000.txt",2)
#print(DiskMass_MW)


# In[5]:


# The mass of component 3 (Bulge stars) for the Milky Way Galaxy is:
#BulgeMass_MW = ComponentMass("MW_000.txt",3)
#print(BulgeMass_MW)


# In[6]:


# Calculating the total mass of the Milky Way Galaxy
#TotalMass_MW = 0 # Setting the initial total mass to 0
# Using a for loop over the particle_type (in the function ComponentMass)
# Prints the total mass of the MW Galaxy
#for i in range(1,4):
 #   TotalMass_MW += ComponentMass("MW_000.txt",i)
#print(TotalMass_MW)


# In[7]:


# Calculating the baryon fraction of MW:
#fbar_MW = (DiskMass_MW+BulgeMass_MW)/TotalMass_MW
#print(np.around(fbar_MW,3))


# In[8]:


# The mass of component 1 (Halo) for Andromeda Galaxy (M31) is:
#HaloMass_M31 = ComponentMass("M31_000.txt",1)
#print(HaloMass_M31)


# In[9]:


# The mass of component 2 (Disk stars) for Andromeda Galaxy (M31) is:
#DiskMass_M31 = ComponentMass("M31_000.txt",2)
#print(DiskMass_M31)


# In[10]:


# The mass of component 3 (Bulge stars) for Andromeda Galaxy (M31) is:
#BulgeMass_M31 = ComponentMass("M31_000.txt",3)
#print(BulgeMass_M31)


# In[11]:


# Calculating the total mass of the Andromeda Galaxy
#TotalMass_M31 = 0 # Setting the initial total mass to 0
# Using a for loop over the particle_type (in the function ComponentMass)
# Prints the total mass of the M31 Galaxy
#for i in range(1,4):
 #   TotalMass_M31 += ComponentMass("M31_000.txt",i)
#print(TotalMass_M31)


# In[12]:


# Calculating the baryon fraction of M31:
#fbar_M31 = (DiskMass_M31+BulgeMass_M31)/TotalMass_M31
#print(np.around(fbar_M31,3))


# In[13]:


# The mass of component 1 (Halo) for the M33 Galaxy is:
#HaloMass_M33= ComponentMass("M33_000.txt",1)
#print(HaloMass_M33)


# In[14]:


# The mass of component 2 (Disk stars) for the M33 Galaxy is:
#DiskMass_M33 = ComponentMass("M33_000.txt",2)
#print(DiskMass_M33)


# In[15]:


# The mass of component 3 (Bulge stars) for the M33 Galaxy is:
#BulgeMass_M33 = ComponentMass("M33_000.txt",3)
#print(BulgeMass_M33)


# In[16]:


# Calculating the total mass of M33 Galaxy
#TotalMass_M33 = 0 # Setting the initial total mass to 0
# Using a for loop over the particle_type (in the function ComponentMass)
# Prints the total mass of the M33 Galaxy
#for i in range(1,4):
 #   TotalMass_M33 += ComponentMass("M33_000.txt",i)
#print(TotalMass_M33)


# In[17]:


# Calculating the baryon fraction of M33:
#fbar_M33 = (DiskMass_M33+BulgeMass_M33)/TotalMass_M33
#print(np.around(fbar_M33,3))


# In[18]:


# Calculating total mass of the local group.
#LocGroupMass = TotalMass_MW + TotalMass_M31 + TotalMass_M33
#print(LocGroupMass)


# In[ ]:




