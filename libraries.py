"""
This file contains useful libraries for BEACON
"""

import sys, os
import numpy as np
import networkx 

def lin2log_metal(Z, eZ):
    FeH = np.log10(Z)-np.log10(0.019)
    eFeH = 0.434*eZ/Z
    return FeH, eFeH

def log2lin_metal(FeH, eFeH):
    Z = 0.019*10**FeH
    eZ = (Z*eFeH)/0.434
    return Z, eZ

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        G.add_nodes_from(part)
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

def rotate_pa(Coo, PA):
    phi = ((90.-PA)/180.)*np.pi

    Rot_matrix = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    Rot_Coo = np.dot(Rot_matrix.T,Coo.T).T

    return Rot_Coo

def deproject_ellip(Coo, Ellipticity = 0., PA = 0.):

   Rot_Coo = rotate_pa(Coo, PA)

   try:
      Rot_Coo[:,1] /= (1-Ellipticity)
   except:
      Rot_Coo[1] /= (1-Ellipticity)

   Dep_Coo = rotate_pa(Rot_Coo, 180.-PA)

   return Dep_Coo
  
def wcs2xy(radec, radec0, Distance_pc):
  
   radec_rad = radec/180.*np.pi
   radec0_rad = radec0/180.*np.pi
 
   xy_ang = np.array([2.*np.arcsin(np.cos(radec_rad[:,1])*np.sin((radec_rad[:,0]-radec0_rad[0])/2.)), radec_rad[:,1]-radec0_rad[1]]).T #HAVERSINE FUNCTION

   xy_pc = Distance_pc*np.tan(xy_ang)
   
   return xy_pc

def xy2wcs(xy_pc, radec0, Distance_pc):
   radec0_rad = radec0/180.*np.pi
   xy_ang = np.arctan(xy_pc/Distance_pc)
   
   try:
      radec_rad = np.array([2.*np.arcsin(np.sin(xy_ang[:,0]/2.)/np.cos(xy_ang[:,1] + radec0_rad[1]))+radec0_rad[0], xy_ang[:,1] + radec0_rad[1]]).T  #INVERTED HAVERSINE FUNCTION
   except:
      radec_rad = np.array([2.*np.arcsin(np.sin(xy_ang[0]/2.)/np.cos(xy_ang[1] + radec0_rad[1]))+radec0_rad[0], xy_ang[1] + radec0_rad[1]]).T  #INVERTED HAVERSINE FUNCTION
   
   radec = radec_rad*180./np.pi

   return radec

def Create_Ellipse(Coo, r, e, PA, Spherical= False, Angle=None):
        
    Coo_RAD = np.zeros(len(Coo))
    Coo_RAD[1] = Coo[1]/180.*np.pi
    Coo_RAD[0] = Coo[0]/180.*np.pi
    
    if Angle is None:
        t = np.arange(0,2.*np.pi,0.05)
        t = np.concatenate((t,[t[0]]),axis=0)
    else:
        t = Angle
       
    phi = ((-PA+90.)/180.)*np.pi
    a = r/180.*np.pi
    b = a*(1.-e)
    
    y_ellipse = Coo_RAD[1]+(a*np.cos(t)*np.sin(phi)+b*np.sin(t)*np.cos(phi))
        
    if Spherical:
        x_non_cor = a*np.cos(t)*np.cos(phi)-b*np.sin(t)*np.sin(phi)
        x_ellipse = 2.*np.arcsin(np.sin(x_non_cor/2.)/np.cos(y_ellipse))+Coo_RAD[0]
    else:
        x_ellipse = Coo_RAD[0] + a*np.cos(t)*np.cos(phi)-b*np.sin(t)*np.sin(phi)
    
    return x_ellipse/np.pi*180., y_ellipse/np.pi*180.

def polar_coo(car_coo):

   pol_coo = np.zeros(car_coo.shape)
   
   pol_coo[:,0] = np.sqrt(car_coo[:,0]**2+car_coo[:,1]**2)
   pol_coo[:,1] = np.arctan2(car_coo[:,1], car_coo[:,0])

   try:
      pol_coo[:,2] = car_coo[:,2]
   except:
     pass
   
   pol_coo[np.where(pol_coo[:,1] < 0), 1] += 2.*np.pi
   
   return pol_coo

def polar_vel(car_coo, car_vel):

   xy2 = car_coo[:,0]**2 + car_coo[:,1]**2
   pol_vel = np.zeros(car_vel.shape)

   pol_vel[:,0] = (car_coo[:,0]/np.sqrt(xy2))*car_vel[:,0] + (car_coo[:,1]/np.sqrt(xy2))*car_vel[:,1]
   pol_vel[:,1] = -(car_coo[:,1]/xy2)*car_vel[:,0] + (car_coo[:,0]/xy2)*car_vel[:,1]
   try:
      pol_vel[:,2] = car_vel[:,2]
   except:
     pass
   
   return pol_vel

def avg_std(values, weights):
    """
    Return the weighted average and standard deviation.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2., weights=weights)
    return [average, np.sqrt(variance)]

def kinetic_pa(XY, V_los, Weights):
    Poss = np.zeros([len(XY),3])
    Poss[:,0:2] = XY
    V = np.zeros([len(XY),3])
    V[:,2] = V_los
    Lp = np.sum(np.cross(Poss, V)*Weights[:, None],axis=0)/np.sum(Weights)
    return Lp

def polar_vector(Angle_Modulus, Galaxy_center):
    Angle, Modulus = Angle_Modulus
   
    x = Galaxy_center[0] + Modulus * np.sin(Angle)
    y = Galaxy_center[1] + Modulus * np.cos(Angle)

    return [x,y]

def Inner_ColorBar(Figure, Axis, Separation_x, Separation_y, Width, Long):
  Cb_coo = [Axis.get_position().get_points()[0,0]+Separation_x, Axis.get_position().get_points()[0,1]+Separation_y, Width, Long]
  Cb_ax =Figure.add_axes(Cb_coo)

  return Cb_ax

"""
Andres del Pino Molina
"""
