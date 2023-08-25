import scipy.io as sio
import numpy as np
import pandas as pd
import datetime, time
import numpy.linalg as lia
import math
import numpy.matlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate


class model():
    
    def __init__(self, setup, model_str, scale):
        self.s = setup
        self.y0 = setup.y0
        self.par = setup.par
        self.i = setup.i
                
        self.t = np.arange(0,  self.s.tf.toordinal() + 1 - self.s.ti.toordinal())

        
        self.dy = np.zeros((6 *len(self.s.geo.popnodes)))
        
        
        fluxes=np.dot(np.exp(-self.s.geo.dist_road/self.p.D), np.diag(self.s.geo.popnodes_ws)) 
        np.fill_diagonal(fluxes,0)
        self.fluxes = fluxes/np.matlib.repmat(np.sum(fluxes,axis=1),365,1)

        if (scale == 'dept'):
            self.fluxes = self.s.geo.ws_adm1.T @ self.fluxes @ self.s.geo.ws_adm1
            np.fill_diagonal(self.fluxes, 0)
            row_sums = self.fluxes.sum(axis=1)
            self.fluxes = self.fluxes / row_sums[:, np.newaxis]
        

        self.model = self.SIB_norm
        
    def run(self):
        tic = time.time()
        sib = scipy.integrate.solve_ivp(fun = self.model,
                                        t_span = (self.t[0], self.t[-1]+1),
                                        y0 = self.y0,
                                        method='RK45',
                                        t_eval = self.t,
                                        max_step = 1.0)
        
        print(">>> Simulation done in ", time.time()-tic)
        return sib
    

    def ODEmodel(self, t, y):
        t = int(t)
        
        dotM = np.dot(self.s.M,  np.dot(self.s.XR, y[self.i.CF]))
        dotE = np.dot(self.WS, y[self.i.ES])
        dotFU = np.dot(self.HU, y[self.i.CF])
        dotFD = np.dot(self.HD, y[self.i.CF])
        
        
        self.dy[self.i.WH] = self.par.beta_FH * dotM - (self.par.muH + self.par.gamma) * y[self.i.WH]
        self.dy[self.i.ES] = self.par.beta_HS * self.par.lambdaE * dotE - (self.par.muS) * y[self.i.ES]
        self.dy[self.i.CF] = self.par.beta_SF * y[self.i.ES] + \
                            self.par.lambdaFU * dotFU + \
                            self.par.lambdaFD * dotFD - \
                            (self.par.muF+self.s.XR+self.par.lambdaFU+self.par.lambdaFD) * y[self.i.CF]
        
        return self.dy 
