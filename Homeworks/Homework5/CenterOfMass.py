#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[2]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities 
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.dot(a,m)/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.dot(b,m)/np.sum(m)
        # zcomponent Center of mass
        Ccom = np.dot(c,m)/np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("maxR", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = self.x - XCOM2 
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [np.around(XCOM,2), np.around(YCOM,2), np.around(ZCOM,2)]

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        return  COMP 
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX
        yV = self.y - COMY
        zV = self.z - COMZ
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV <= RVMAX.value)

        # determine the velocity and mass of those particles within the max radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [np.around(VXCOM,2), np.around(VYCOM,2), np.around(VZCOM,2)]

        # return the COM vector                                                                                        
        return COMV 


# In[3]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)


# In[4]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COMP = MWCOM.COM_P(0.1) 
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2]) 
print(MW_COMP *u.kpc)
print(MW_COMV *u.km/u.s)


# In[5]:


M31COM = CenterOfMass("M31_000.txt", 2)

#M31:store the position and velocity COM
M31_COMP =M31COM.COM_P(0.1) 
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
print(M31_COMP *u.kpc)
print(M31_COMV *u.km/u.s)


# In[6]:


M33COM = CenterOfMass("M33_000.txt", 2)

#M33: store the position and velocity COM
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])
print(M33_COMP *u.kpc)
print(M33_COMV *u.km/u.s)


# In[7]:


# now write your own code to answer questions
# The length of the center of mass vector of Milky Way:
COM_MW = np.around(np.linalg.norm(MW_COMP),2) * u.kpc
# The length of the center of mass vector of M31:
COM_M31 = np.around(np.linalg.norm(M31_COMP),2) * u.kpc
# The distance between MW and M31 is: 
print(COM_M31 - COM_MW)


# In[8]:


# The speed of the MW is:
MW_vel = np.around(np.linalg.norm(MW_COMV),2) * u.km/u.s
# The speed of the Andromeda galaxy is: 
M31_vel = np.around(np.linalg.norm(M31_COMV),2) * u.km/u.s
# The relative speed of M31 with respect to MW is:
print (np.absolute((M31_vel - MW_vel)))


# In[9]:


# the seperation between M31 and M33
sep_M31_M33 = np.array(M33_COMP) - np.array(M31_COMP) 
# the magnitude of the distance is
print(np.around(np.linalg.norm(sep_M31_M33),2) * u.kpc) 


# In[10]:


# the relative velocity between M31 and M33
vel_M31_M33 = np.array(M33_COMV) - np.array(M31_COMV) 
# the magnitude of the velocity is:
print(np.around(np.linalg.norm(vel_M31_M33),2)*u.km/u.s)


# The determination of the COM of the MW and M31 is an iterative process because we need to know if they will experience any significant deviations (due to influences from other sources) from their collision course in the coming years. This can be done only by measuring the COM of MW and M31 at multiple (astronomical) timesteps. 

# In[ ]:




