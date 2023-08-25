#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:13:05 2023

@author: cristiano
"""

import numpy as np
import pandas as pd
import scipy.io as sio

class OPIsetup:
    
    def model(self):
        
        datapath = '../dataOCN/'

        self.dt = 1

        self.par = Parameters.from_list([1.1557,        #theta  CALI
                                            0.1772,        #m  CALI
                                            98.5159,       #D  CALI
                                            0.4663,        #phi CALI
                                            0.0523,        #sigma
                                            0.0265,        #mub 
                                            0.2344,        #beta01
                                            0.2344,        #beta02
                                            0.2344,        #beta03
                                            0.2344,        #beta04
                                            0.2344,        #beta05
                                            0.2344,        #beta06
                                            0.2344,        #beta07
                                            0.2344,        #beta08
                                            0.2344,        #beta09
                                            0.2344,        #beta010
                                            0.5])           #r        
        
        # Read geography data
        self.geo = Geography(datapath)
        
        
        # Setup
        self.i = OPI_Indexing(self.geo.nNodes)
        self.y0 = self.Build_Initial_Condtion()
        
        self.geo.nNodes_ws = self.geo.nNodes
        self.geo.popnodes_ws =self.geo.popnodes
        
                
    def Build_Initial_Condtion(self):
        ic_pts = [150-1, 109-1, 119-1, 118-1]
        
        I0 = np.round(np.minimum(
            np.insert(1100*self.geo.popnodes.take(ic_pts[1:])/np.sum(self.geo.popnodes.take(ic_pts[1:])),0,1000)\
                    /(self.p.sigma*self.p.rep_ratio),
                self.geo.popnodes.take(ic_pts)))   # Tested OK
        
        y0 = np.zeros((6 *len(self.geo.popnodes)))
        y0[self.i.S]= self.geo.popnodes.copy()
        np.put(y0, ic_pts, y0[self.i.S].take(ic_pts) - I0)             # Tested OK
        np.put(y0, ic_pts + self.i.I[0], np.around(self.p.sigma*I0))                   # Tested OK
        np.put(y0, ic_pts + self.i.A[0], np.around((1-self.p.sigma)*I0))              # Tested OK
        np.put(y0, ic_pts + self.i.R[0], 0)
        np.put(y0, ic_pts + self.i.B[0], (self.p.r * y0[self.i.A].take(ic_pts) + y0[self.i.I].take(ic_pts))*self.p.theta/self.geo.popnodes.take(ic_pts)/self.p.muB)   # Tested OK
        np.put(y0, ic_pts + self.i.C[0], np.around(self.p.sigma*I0))                   # Tested OK
        
        return y0


class Results():
    def __init__(self, ODEmodel, setup):
        self.I = pd.DataFrame(ODEmodel[setup.i.I,:].T,columns=np.arange(0,setup.geo.nNodes), index=pd.date_range(setup.t1i, setup.t2f))
        self.S = pd.DataFrame(ODEmodel[setup.i.S,:].T,columns=np.arange(0,setup.geo.nNodes), index=pd.date_range(setup.t1i, setup.t2f))
        self.F = pd.DataFrame(ODEmodel[setup.i.F,:].T,columns=np.arange(0,setup.geo.nNodes), index=pd.date_range(setup.t1i, setup.t2f))



class Parameters():
    """ VERY IMPORTANT!!!!!
    Before continuing, remember to modify the param_list vector"""
    def __init__(self, theta, m, D, phi, sigma, muB, beta0, r):
        self.theta  = theta
        self.m      = m
        self.D      = D
        self.phi    = phi
        self.sigma  = sigma
        self.muB    = muB
        self.beta0  = beta0
        self.r      = r
        
        
        # Fixed Parameters
        self.gamma     = 0.2                 # rate at which people recover from cholera (day^-1)
        self.rep_ratio = 1                   # reported to total (symptomatic) cases ratio --> asymptomatic cases do not get reported

        # Parameters fixed according to litterature and former research
        self.mu          = 1/(61.4*365)              # population natality and mortality rate (day^-1)
        self.alpha       = -np.log(0.98)*self.gamma  # mortality rate due to cholera (day^-1) (2% case fatality rate)
        self.epsilon_muB = 0                         # possible seasonality in mortality of bacteria (set to 0)
        
         #From code SERRA
        self.muB_wsd = np.full((10, ), 0.2) 
    @classmethod    
    def from_list(cls, param_list):
        """ Maintains compatility with enrico's prior"""
        return cls(    param_list[0],
                       param_list[1],
                       param_list[2],
                       param_list[3],
                       param_list[4],
                       param_list[5],
                       param_list[6:16],
                       param_list[16])
        
class Geography():
    def __init__(self, directory):
        subcatchments = sio.loadmat(directory+'subcatchments.mat')


        CTC = subcatchments['CTC'].astype('int')
        X = subcatchments['X']
        Y = subcatchments['Y']

        NX = np.ceil(np.max(X)).astype('int')
        NY = np.ceil(np.max(Y)).astype('int')
        
        FD = sio.loadmat(directory+'FD.mat')
        
        self.FD_A = FD['FD_A'].flatten()
        self.FD_X = FD['FD_X'].flatten()
        self.FD_Y = FD['FD_Y'].flatten()
        self.FD_downNode = FD['FD_downNode'].astype('int').flatten()

        self.FD_nNodes = len(self.FD_A)
        
        
        

        self.nNodes = geodata['nNodes'].flatten()[0]
        self.x = geodata['X']
        self.y = geodata['Y']
        self.AD = geodata['AD']
        self.popnodes = geodata['POPnodes'].flatten()
        self.ws_adm1 = geodata['WS_dept']
        self.dist_road = geodata['dist_road']
        self.H = np.sum(self.popnodes)
        
    
class OPI_Indexing():
    def __init__ (self, nNodes):
        self.I = np.arange(0*nNodes, 1*nNodes)
        self.S = np.arange(1*nNodes, 2*nNodes)
        self.F = np.arange(2*nNodes, 3*nNodes)