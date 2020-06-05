#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:18:44 2020

@author: Rafael Veiga rafaelvalenteveiga@gmail.com
@author: matheustorquato matheusft@gmail.com
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as spi
import pyswarms as ps
from pyswarms.utils.plotters import plot_cost_history
import pickle as pk

class Models:
    def __init__(self):
        self.isFit=False
        self.isBetaChange = False
        
    def __genBoot(self,series, times = 500):
        series = np.diff(series)
        series = np.insert(series, 0, 1)
        series[series < 0] = 0
        results = []
        for i in range(0,times):
            results.append(np.random.multinomial(n = sum(series), pvals = series/sum(series)))
        return np.array(results)
    
    def __getConfidenceInterval(self, series, length, level):
        series = np.array(series)
    
        #Compute mean value
        meanValue = [np.mean(series[:,i]) for i in range(0,length)]

        #Compute deltaStar
        deltaStar = meanValue - series
        #Compute lower and uper bound
        q= (1-level)/2
        deltaL = [np.quantile(deltaStar[:,i], q = q) for i in range(0,length)]
        deltaU = [np.quantile(deltaStar[:,i], q = 1-q) for i in range(0,length)]

        #Compute CI
        lowerBound  = np.array(meanValue) + np.array(deltaL)
        UpperBound  = np.array(meanValue) + np.array(deltaU)
        return [meanValue, lowerBound, UpperBound]
        
    def getResiduosQuadatico(self):
        y = np.array(self.y)
        d = np.array(self.d)
        ypred = np.array(self.ypred)
        dpred = np.array(self.dpred)
        y = y[0:len(self.x)]
        d = d[0:len(self.x)]
        ypred = ypred[0:len(self.x)]
        dpred = dpred[0:len(self.x)]
        return ((y - ypred)**2)*(1-self.pesoMorte) + ((d-dpred)**2)*self.pesoMorte

    def getReQuadPadronizado(self):
        y = np.array(self.y)
        d = np.array(self.d)
        ypred = np.array(self.ypred)
        dpred = np.array(self.dpred)
        y = y[0:len(self.x)]
        d = d[0:len(self.x)]
        ypred = ypred[0:len(self.x)]
        dpred = dpred[0:len(self.x)]
        return (((y - ypred)**2)/np.sqrt(ypred+1))*(1-self.pesoMorte) + (((d-dpred)**2)/np.sqrt(dpred+1))*self.pesoMorte
    
    def plotCost(self):
        if self.isFit:
            plot_cost_history(cost_history=self.cost_history)
            plt.show()
        else:
            print('\nModels is not fitted\n')

    def save(self,fileName):
        file = open(fileName,'wb')
        pk.dump(self,file)
        file.close()
        
    def load(fileName):
        file = open(fileName,'rb')
        model = pk.load(file)
        return model

class SIR:
    ''' SIR Model'''
    def __init__(self,tamanhoPop,numeroProcessadores=None):
        super(SIR,self).__init__()
        self.N = tamanhoPop
        self.numeroProcessadores = numeroProcessadores
    
    def __cal_EDO(self,x,beta,gamma):
            ND = len(x)-1
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            beta = np.array(beta)
            gamma = np.array(gamma)
            def SIR_diff_eqs(INP, t, beta, gamma):
                Y = np.zeros((3))
                V = INP
                Y[0] = - beta * V[0] * V[1]                 #S
                Y[1] = beta * V[0] * V[1] - gamma * V[1]    #I
                Y[2] = gamma * V[1]                         #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (self.S0, self.I0,self.R0), t_range,
                                    args=(beta, gamma))
            
            S=result_fit[:, 0]*self.N
            R=result_fit[:, 2]*self.N
            I=result_fit[:, 1]*self.N
            
            return S,I,R
        
    def __cal_EDO_2(self,x,beta1,gamma,beta2,tempo):
            ND = len(x)-1
            
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            def H(t):
                h = 1.0/(1.0+ np.exp(-2.0*50*t))
                return h
            def beta(t,t1,b,b1):
                beta = b*H(t1-t) + b1*H(t-t1) 
                return beta

            gamma = np.array(gamma)
            def SIR_diff_eqs(INP, t, beta1, gamma,beta2,t1):
                Y = np.zeros((3))
                V = INP
                Y[0] = - beta(t,t1,beta1,beta2) * V[0] * V[1]                 #S
                Y[1] = beta(t,t1,beta1,beta2) * V[0] * V[1] - gamma * V[1]    #I
                Y[2] = gamma * V[1]                         #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (self.S0, self.I0,self.R0), t_range,
                                    args=(beta1, gamma,beta2,tempo))
            
            S=result_fit[:, 0]*self.N
            R=result_fit[:, 2]*self.N
            I=result_fit[:, 1]*self.N
            
            return S,I,R
    
    def objectiveFunction(self,coef,x ,y,stand_error):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        y = y*self.N
        if stand_error:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3])
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
            elif self.isBetaChange:
                for i in range(tam2):
                    S,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange)
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
            else:
                for i in range(tam2):
                    S,I,R = self.__cal_EDO(x,coef[i,0],coef[i,1])
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
        else:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3])
                    soma[i]= (((y-(I+R)))**2).mean()
            elif self.isBetaChange:
                for i in range(tam2):
                    S,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange)
                    soma[i]= (((y-(I+R)))**2).mean()
            else:
                for i in range(tam2):
                    S,I,R = self.__cal_EDO(x,coef[i,0],coef[i,1])
                    soma[i]= (((y-(I+R)))**2).mean()
        return soma
    def fit(self, y , bound = ([0,1/21],[1,1/5]),stand_error=True, isBetaChange=False, dayBetaChange = None,particles=50,itera=500,c1= 0.5, c2= 0.3, w = 0.9, k=3,p=1):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        x = range(1,len(y)+1)
        self.isBetaChange = isBetaChange
        self.dayBetaChange = dayBetaChange
        self.y = y
        self.x = x
        df = np.array(y)/self.N
        self.I0 = df[0]
        self.S0 = 1-self.I0
        self.R0 = 0
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':p}
        optimizer = None
        if bound==None:
            if (isBetaChange) & (dayBetaChange==None):
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=4, options=options)
            elif isBetaChange:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=3, options=options)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=2, options=options)                
        else:
            if (isBetaChange) & (dayBetaChange==None):
                if len(bound[0])==2:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].append(bound[0][0])
                    bound[1].append(bound[1][0])
                    bound[0].append(x[4])
                    bound[1].append(x[-5])
                    bound[0][3] = x[4]
                    bound[1][3] = x[-5]
                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=4, options=options,bounds=bound)
            elif isBetaChange:
                if len(bound[0])==2:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].append(bound[0][1])
                    bound[1].append(bound[1][1])
                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=3, options=options,bounds=bound)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=2, options=options,bounds=bound)
                
        cost = pos = None
        if isBetaChange:
            cost, pos = optimizer.optimize(self.objectiveFunction, itera, x = x,y=df,stand_error=stand_error,n_processes=self.numeroProcessadores)
        else:
            cost, pos = optimizer.optimize(self.objectiveFunction, itera, x = x,y=df,stand_error=stand_error,n_processes=self.numeroProcessadores)
            self.beta = pos[0]
            self.gamma = pos[1]
        if isBetaChange:
            self.beta1 = pos[0]
            self.gamma = pos[1]
            self.beta2 = pos[2]
            if dayBetaChange==None:
                self.dayBetaChange = pos[3]
            else:
                self.dayBetaChange = dayBetaChange
        self.rmse = cost
        self.cost_history = optimizer.cost_history
            
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays) 
        
        if self.isBetaChange:
            S,I,R = self.__cal_EDO_2(x,self.beta1,self.gamma,self.beta2,self.dayBetaChange)
        else:
            S,I,R = self.__cal_EDO(x,self.beta,self.gamma)
        self.ypred = I+R
        self.S = S
        self.I = I
        self.R = R         
        return self.ypred

    def plot(self,local):
        ypred = self.predict(self.x)
        plt.plot(ypred,c='b',label='Predição Infectados')
        plt.plot(self.y,c='r',marker='o', markersize=3,label='Infectados')
        plt.legend(fontsize=15)
        plt.title('Dinâmica do CoviD19 - {}'.format(local),fontsize=20)
        plt.ylabel('Casos COnfirmados',fontsize=15)
        plt.xlabel('Dias',fontsize=15)
        plt.show()
    def getCoef(self):
        if self.isBetaChange:
            return ['beta1','beta2','gamma','dayBetaChange'],[self.beta1,self.beta2,self.gamma,self.dayBetaChange]
        return ['beta','gamma'], [self.beta,self.gamma]

    def plotFit(self):
        plt.style.use('seaborn-deep')
        fig, axes = plt.subplots(figsize = (18,8))
        try:
            plt.plot(self.x, self.ypred, label = "Fitted", c = "red")
            plt.scatter(self.x, self.y, label = "Observed", c = "blue")
            plt.legend(loc='upper left')
            plt.show()
        except:
            print("There is no predicted value")

class SEIR:
    ''' SIR Model'''
    def __init__(self,tamanhoPop,numeroProcessadores=None):
        super(SEIR,self).__init__()
        self.N = tamanhoPop
        self.numeroProcessadores = numeroProcessadores
    
    def __cal_EDO(self,x,beta,gamma,mu,sigma):
            ND = len(x)-1
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            #beta = np.array(beta)
            #gamma = np.array(gamma)
            #mu = np.array(mu)
            #sigma = np.array(sigma)
            
            def SEIR_diff_eqs(INP, t, beta, gamma,mu,sigma):
                Y = np.zeros((4))
                V = INP
                Y[0] = mu - beta * V[0] * V[2] - mu * V[0]  # Susceptile
                Y[1] = beta * V[0] * V[2] - sigma * V[1] - mu * V[1] # Exposed
                Y[2] = sigma * V[1] - gamma * V[2] - mu * V[2] # Infectious
                Y[3] = gamma * V[2] #recuperado
                return Y   # For odeint

                return Y
            result_fit = spi.odeint(SEIR_diff_eqs, (self.S0,self.E0, self.I0,self.R0), t_range,
                                    args=(beta, gamma,mu,sigma))
            
            S=result_fit[:, 0]*self.N
            E=result_fit[:, 1]*self.N
            I=result_fit[:, 2]*self.N
            R=result_fit[:, 3]*self.N
            
            return S,E,I,R
      
    def __cal_EDO_2(self,x,beta1,beta2,dayBetaChange,gamma,mu,sigma):
            ND = len(x)-1
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            #beta1 = np.array(beta1)
            #beta2 = np.array(beta2)
            #gamma = np.array(gamma)
            #mu = np.array(mu)
            #sigma = np.array(sigma)
            def Hf(t):
                h = 1.0/(1.0+ np.exp(-2.0*50*t))
                return h
            def beta(t,t1,b,b1):
                beta = b*Hf(t1-t) + b1*Hf(t-t1) 
                return beta
            def SEIR_diff_eqs(INP, t, beta1,beta2,t1, gamma,mu,sigma):
                Y = np.zeros((4))
                V = INP
                Y[0] = mu - beta(t,t1,beta1,beta2) * V[0] * V[2] - mu * V[0]  # Susceptile
                Y[1] = beta(t,t1,beta1,beta2) * V[0] * V[2] - sigma * V[1] - mu * V[1] # Exposed
                Y[2] = sigma * V[1] - gamma * V[2] - mu * V[2] # Infectious
                Y[3] = gamma * V[2] #recuperado
                return Y   # For odeint

                return Y
            result_fit = spi.odeint(SEIR_diff_eqs, (self.S0,self.E0, self.I0,self.R0), t_range,
                                    args=(beta1,beta2,dayBetaChange, gamma,mu,sigma))
            
            S=result_fit[:, 0]*self.N
            E=result_fit[:, 1]*self.N
            I=result_fit[:, 2]*self.N
            R=result_fit[:, 3]*self.N
            
            return S,E,I,R
    def __objectiveFunction(self,coef,x ,y,stand_error):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        #__cal_EDO(self,x,beta,gamma,mu,sigma)
        #__cal_EDO2(self,x,beta1,beta2,dayBetaChange,gamma,mu,sigma)
        if stand_error:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],self.mu,coef[i,4])
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],self.dayBetaChange,coef[i,2],self.mu,coef[i,3])
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
            else:
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO(x,coef[i,0],coef[i,1],self.mu,coef[i,2])
                    soma[i]= (((y-(I+R))/np.sqrt((I+R)+1))**2).mean()
        else:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],self.mu,coef[i,4])
                    soma[i]= (((y-(I+R)))**2).mean()
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO_2(x,coef[i,0],coef[i,1],self.dayBetaChange,coef[i,2],self.mu,coef[i,3])
                    soma[i]= (((y-(I+R)))**2).mean()
            else:
                for i in range(tam2):
                    S,E,I,R = self.__cal_EDO(x,coef[i,0],coef[i,1],self.mu,coef[i,2])
                    soma[i]= (((y-(I+R)))**2).mean()
        return soma
    

    def fit(self, y , bound = ([0,1/7,1/6],[1.5,1/4,1/4]) ,stand_error=True, isBetaChange=True,dayBetaChange = None,particles=50,itera=500,c1=0.3,c2= 0.3, w= 0.9,k=3,p=2):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        x = range(1,len(y)+1)
        self.y = y
        self.I0 = np.array(y[0])/self.N
        self.S0 = 1-self.I0
        self.R0 = 0
        self.E0 = 0
        self.mu = 1/(75.51*365)
        # options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
        # if bound==None:
        #     optimizer = ps.single.GeneralOptimizerPSO(n_particles=50, dimensions=3, options=options,topology=Star())
        #     cost, pos = optimizer.optimize(self.__objectiveFunction, 500, x = x,y=y,mu=1/(75.51*365),n_processes=self.numeroProcessadores)
        #     self.beta = pos[0]
        #     self.gamma = pos[1]
        #     self.mu = 1/(75.51*365)
        #     self.sigma = pos[2]
        #     self.x = x
        #     self.rmse = cost
        #     self.optimize = optimizer
            
        # else:
        #     optimizer = ps.single.GeneralOptimizerPSO(n_particles=50, dimensions=3, options=options,bounds=bound,topology=Star())
        #     cost, pos = optimizer.optimize(self.__objectiveFunction, 500, x = x,y=y,mu=1/(75.51*365),n_processes=self.numeroProcessadores)
        #     self.beta = pos[0]
        #     self.gamma = pos[1]
        #     self.mu = 1/(75.51*365)
        #     self.sigma = pos[2]
        #     self.x = x
        #     self.rmse = cost
        #     self.optimize = optimizer
        self.isBetaChange = isBetaChange
        self.dayBetaChange = dayBetaChange
        self.y = y
        self.x = x

        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':p}
        optimizer = None
        
        if bound==None:
            if (isBetaChange) & (dayBetaChange==None):
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=5, options=options)
            elif isBetaChange:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=4, options=options)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=3, options=options)                
        else:
            if (isBetaChange) & (dayBetaChange==None):
                if len(bound[0])==3:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].insert(1,bound[0][0])
                    bound[1].insert(1,bound[1][0])
                    bound[0].insert(2,x[4])
                    bound[1].insert(2,x[-5])

                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=5, options=options,bounds=bound)
            elif isBetaChange:
                if len(bound[0])==3:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].insert(1,bound[0][0])
                    bound[1].insert(1,bound[1][0])
                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=4, options=options,bounds=bound)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=3, options=options,bounds=bound)
                
        cost = pos = None
        if isBetaChange:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=y,stand_error=stand_error,n_processes=self.numeroProcessadores)
        else:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=y,stand_error=stand_error,n_processes=self.numeroProcessadores)
            self.beta = pos[0]
            self.gamma = pos[1]
            self.sigma = pos[2]
            
        if isBetaChange:
            self.beta1 = pos[0]
            self.beta2 = pos[1]
            
            if dayBetaChange==None:
                self.dayBetaChange = pos[2]
                self.gamma = pos[3]
                self.sigma = pos[4]
            else:
                self.dayBetaChange = dayBetaChange
                self.gamma = pos[2]
                self.sigma = pos[3]

        self.rmse = cost
        self.cost_history = optimizer.cost_history
            
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays)
        
        if self.isBetaChange:
            S,E,I,R = self.__cal_EDO_2(x,self.beta1,self.beta2,self.dayBetaChange,self.gamma,self.mu,self.sigma)
        else:
            S,E,I,R = self.__cal_EDO(x,self.beta,self.gamma,self.mu,self.sigma)
        self.ypred = I+R
        self.S = S
        self.E = E
        self.I = I
        self.R = R         
        return self.ypred

    def plot(self,local):
        ypred = self.predict(self.x)
        plt.plot(ypred,c='b',label='Predição Infectados')
        plt.plot(self.y,c='r',marker='o', markersize=3,label='Infectados')
        plt.legend(fontsize=15)
        plt.title('Dinâmica do CoviD19 - {}'.format(local),fontsize=20)
        plt.ylabel('Casos COnfirmados',fontsize=15)
        plt.xlabel('Dias',fontsize=15)
        plt.show()
    def getCoef(self):
        #__cal_EDO(self,x,beta,gamma,mu,sigma)
        #__cal_EDO2(self,x,beta1,beta2,dayBetaChange,gamma,mu,sigma)
        if self.isBetaChange:
            return ['beta1','beta2','dayBetaChange','gamma','mu','sigma'],[self.beta1,self.beta2,self.dayBetaChange,self.gamma,self.mu,self.sigma]
        return ['beta','gamma','mu','sigma'],[self.beta,self.gamma,self.mu,self.sigma]

    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        #Define empty lists to recive results
        self.__lypred = []
        self.__lspred = []
        self.__lepred = []
        self.__lipred = []
        self.__lrpred = []
        if self.isBetaChange:
            self.__lbeta1 = []
            self.__lbeta2 = []
        else:
            self.__lbeta=[]
        
        self.__lgammaH = []
        self.__lgammaU = []
        self.__ldelta = []
        self.__le0 = []
        self.__lia0 = []
        self.__lis0 = []

        self.y = self.y[np.argsort(self.y)]
        self.d = self.d[np.argsort(self.d)]
        casesSeries = self.__genBoot(self.y, times)
        deathSeries = self.__genBoot(self.d, times)
        for i in range(0,len(casesSeries)):
            super().fit(x = range(1,len(self.y) + 1),
                        y = casesSeries[i],
                        d = deathSeries[i])
            super().predict(self.x)

            self.__lypred.append(self.ypred)
            self.__ldpred.append(self.dpred)
            self.__lhpred.append(self.H)
            self.__lupred.append(self.U)
            self.__lspred.append(self.S)
            self.__lepred.append(self.E)
            self.__lrpred.append(self.R)
            self.__lIA.append(self.IA)
            self.__lIS.append(self.IS)

            self.__lbeta1.append(self.beta1)
            self.__lbeta2.append(self.beta2)
            self.__lgammaH.append(self.gammaH)
            self.__lgammaU.append(self.gammaU)
            self.__ldelta.append(self.delta)
            self.__le0.append(self.e0)
            self.__lia0.append(self.ia0)
            self.__lis0.append(self.is0)
            self.ciResults = {}
        for i, serie in zip([self.__lypred, self.__ldpred, self.__lhpred, self.__lupred, self.__lspred, 
                  self.__lepred, self.__lrpred, self.__lIA, self.__lIS],
                  ["ypred", "dpred", "hpred", "upred", "spred",
                   "epred", "dpred", "IA", "IS"]):
            meanValue, lb, ub = self.__getConfidenceInterval(i,len(i[0]))
            self.ciResults.update({serie: {"meanValue":meanValue,
                 "lb": lb,
                 "ub": ub}})


        #return df with parameters
        self.parameters = pd.DataFrame.from_dict({
            "beta1": self.__lbeta1,
            "beta2": self.__lbeta2,
            "gammaH": self.__lgammaH,
            "gammaU": self.__lgammaU,
            "delta": self.__ldelta,
            "e0": self.__le0,
            "a0": self.__lia0,
            "is0": self.__lis0
        })
    
class SEIRHUD(Models):
    ''' SEIRHU Model'''
    def __init__(self,tamanhoPop,numeroProcessadores=None):
        super(SEIRHUD,self).__init__()
        self.N = tamanhoPop
        self.numeroProcessadores = numeroProcessadores
    
    def __cal_EDO(self,x,beta,gammaH,gammaU,delta,h,ia0,is0,e0):
            ND = len(x)-1
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            beta = np.array(beta)
            delta = np.array(delta)
            def SIR_diff_eqs(INP, t, beta,gammaH,gammaU, delta,h):
                Y = np.zeros((9))
                V = INP
                Y[0] = - beta*V[0]*(V[3] + delta*V[2])                    #S
                Y[1] = beta*V[0]*(V[3] + delta*V[2]) -self.kappa * V[1]
                Y[2] = (1-self.p)*self.kappa*V[1] - self.gammaA*V[2]
                Y[3] = self.p*self.kappa*V[1] - self.gammaS*V[3]
                Y[4] = h*self.xi*self.gammaS*V[3] + (1-self.muU + self.omegaU*self.muU)*gammaU*V[5] -gammaH*V[4]
                Y[5] = h*(1-self.xi)*self.gammaS*V[3] +self.omegaH*gammaH*V[4] -gammaU*V[5]
                Y[6] = self.gammaA*V[2] + (1-(self.muH))*(1-self.omegaH)*gammaH*V[4] + (1-h)*self.gammaS*V[3]
                Y[7] = (1-self.omegaH)*self.muH*gammaH*V[4] + (1-self.omegaU)*self.muU*gammaU*V[5]#R
                Y[8] = self.p*self.kappa*V[1] 
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (1-ia0-is0-e0,e0 ,ia0,is0,0,0,0,0,0), t_range,
                                    args=(beta,gammaH,gammaU, delta,h))
            
            S=result_fit[:, 0]*self.N
            E = result_fit[:, 1]*self.N
            IA=result_fit[:, 2]*self.N
            IS=result_fit[:, 3]*self.N
            H=result_fit[:, 4]*self.N
            U=result_fit[:, 5]*self.N
            R=result_fit[:, 6]*self.N
            D=result_fit[:, 7]*self.N
            Nw=result_fit[:, 8]*self.N
            
            return S,E,IA,IS,H,U,R,D,Nw
        
    def __cal_EDO_2(self,x,beta1,beta2,tempo,gammaH,gammaU,delta,h,ia0,is0,e0):
            ND = len(x)-1
            
            t_start = 0.0
            t_end = ND
            t_inc = 1
            t_range = np.arange(t_start, t_end + t_inc, t_inc)
            def Hf(t):
                h = 1.0/(1.0+ np.exp(-2.0*50*t))
                return h
            def beta(t,t1,b,b1):
                beta = b*Hf(t1-t) + b1*Hf(t-t1) 
                return beta

            delta = np.array(delta)
            def SIR_diff_eqs(INP, t, beta1, beta2,t1,gammaH,gammaU, delta,h):
                #Y[0] = - beta(t,t1,beta1,beta2) * V[0] * V[1]                 #S
                Y = np.zeros((9))
                V = INP
                Y[0] = - beta(t,t1,beta1,beta2)*V[0]*(V[3] + delta*V[2])                    #S
                Y[1] = beta(t,t1,beta1,beta2)*V[0]*(V[3] + delta*V[2]) -self.kappa * V[1]
                Y[2] = (1-self.p)*self.kappa*V[1] - self.gammaA*V[2]
                Y[3] = self.p*self.kappa*V[1] - self.gammaS*V[3]
                Y[4] = h*self.xi*self.gammaS*V[3] + (1-self.muU + self.omegaU*self.muU)*gammaU*V[5] -gammaH*V[4]
                Y[5] = h*(1-self.xi)*self.gammaS*V[3] +self.omegaH*gammaH*V[4] -gammaU*V[5]
                Y[6] = self.gammaA*V[2] + (1-(self.muH))*(1-self.omegaH)*gammaH*V[4] + (1-h)*self.gammaS*V[3]
                Y[7] = (1-self.omegaH)*self.muH*gammaH*V[4] + (1-self.omegaU)*self.muU*gammaU*V[5]#R
                Y[8] = self.p*self.kappa*V[1]                      #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (1-ia0-is0-e0,e0 ,ia0,is0,0,0,0,0,0), t_range,
                                    args=(beta1,beta2,tempo,gammaH,gammaU, delta,h))
            
            S=result_fit[:, 0]*self.N
            E = result_fit[:, 1]*self.N
            IA=result_fit[:, 2]*self.N
            IS=result_fit[:, 3]*self.N
            H=result_fit[:, 4]*self.N
            U=result_fit[:, 5]*self.N
            R=result_fit[:, 6]*self.N
            D=result_fit[:, 7]*self.N
            Nw=result_fit[:, 8]*self.N
            
            return S,E,IA,IS,H,U,R,D,Nw
    
    def objectiveFunction(self,coef,x ,y,d,stand_error):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        if stand_error:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8],coef[i,9])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*(1-self.pesoMorte)+(((d-(D))/np.sqrt(D+1))**2).mean()*self.pesoMorte
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange,coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*(1-self.pesoMorte)+(((d-(D))/np.sqrt(D+1))**2).mean()*self.pesoMorte
            else:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*(1-self.pesoMorte)+(((d-(D))/np.sqrt(D+1))**2).mean()*self.pesoMorte
        else:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8],coef[i,9])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange,coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
            else:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
        return soma
    def fit(self, y, d, pesoMorte = 0.5, kappa = 1/4,p = 0.2,gammaA = 1/3.5, gammaS = 1/4.001, muH = 0.15,
            muU = 0.4,xi = 0.53,omegaU = 0.29,omegaH=0.14 , bound = [[0,1/8,1/12,0,0],[2,1/4,1/3,0.7,0.35]],
            stand_error = True, isBetaChange = False, dayBetaChange = None, particles = 300, itera = 1000, c1 = 0.1, c2 = 0.3, w = 0.9, k = 5, norm = 2):
        '''
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        if len(y)!=len(d):
            print('\ny and d must have the same length\n')
            return
        x = range(1,len(y)+1)
        
        if len(bound)==2:
            if len(bound[0])==5:
                bound[0]=bound[0].copy()
                bound[1]=bound[1].copy()
                bound[0].append(0)
                bound[0].append(0)
                bound[0].append(0)
                bound[1].append(10/self.N)
                bound[1].append(10/self.N)
                bound[1].append(10/self.N)
        self.pesoMorte = pesoMorte
        self.kappa = kappa
        self.p = p
        self.gammaA = gammaA
        self.gammaS = gammaS
        self.muH = muH
        self.muU = muU
        self.xi = xi
        self.omegaU = omegaU
        self.omegaH = omegaH
        self.isBetaChange = isBetaChange
        self.dayBetaChange = dayBetaChange
        self.y = y
        self.d = d
        self.x = x
        df = np.array(y)
        dd = np.array(d)
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
        optimizer = None
        if bound==None:
            if (isBetaChange) & (dayBetaChange==None):
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=10, options=options)
            elif isBetaChange:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=9, options=options)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=8, options=options)                
        else:
            if (isBetaChange) & (dayBetaChange==None):
                if len(bound[0])==8:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].insert(1,bound[0][0])
                    bound[1].insert(1,bound[1][0])
                    bound[0].insert(2,x[4])
                    bound[1].insert(2,x[-5])

                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=10, options=options,bounds=bound)
            elif isBetaChange:
                if len(bound[0])==8:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].insert(1,bound[0][0])
                    bound[1].insert(1,bound[1][0])
                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=9, options=options,bounds=bound)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=8, options=options,bounds=bound)
                
        cost = pos = None
        #__cal_EDO(self,x,beta,gammaH,gammaU,delta,h,ia0,is0,e0)
        #__cal_EDO_2(self,x,beta1,beta2,tempo,gammaH,gammaU,delta,h,ia0,is0,e0)
        if isBetaChange:
            #cost, pos = optimizer.optimize(self.objectiveFunction,itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores, verbose = True)
            cost, pos = optimizer.optimize(self.objectiveFunction,itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores)
        else:
            #cost, pos = optimizer.optimize(self.objectiveFunction, itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores, verbose = True)
            cost, pos = optimizer.optimize(self.objectiveFunction, itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores)
            self.beta = pos[0]
            self.gammaH = pos[1]
            self.gammaU = pos[2]
            self.delta = pos[3]
            self.h = pos[4]
            self.ia0 = pos[5]
            self.is0 = pos[6]
            self.e0 = pos[7]
        if isBetaChange:
            self.beta1 = pos[0]
            self.beta2 = pos[1]
            
            if dayBetaChange==None:
                self.dayBetaChange = pos[2]
                self.gammaH = pos[3]
                self.gammaU = pos[4]
                self.delta = pos[5]
                self.h = pos[6]
                self.ia0 = pos[7]
                self.is0 = pos[8]
                self.e0 = pos[9]
            else:
                self.dayBetaChange = dayBetaChange
                self.gammaH = pos[2]
                self.gammaU = pos[3]
                self.delta = pos[4]
                self.h = pos[5]
                self.ia0 = pos[6]
                self.is0 = pos[7]
                self.e0 = pos[8]
        self.rmse = cost
        self.cost_history = optimizer.cost_history
        self.isFit=True

    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays)   
        self.x = x
        if self.isBetaChange:
            S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,self.beta1,self.beta2,self.dayBetaChange,self.gammaH,self.gammaU,self.delta,self.h,self.ia0,self.is0,self.e0)
        else:
            S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,self.beta,self.gammaH,self.gammaU,self.delta,self.h,self.ia0,self.is0,self.e0)
        self.ypred = Nw
        self.dpred = D
        self.S = S
        self.E = E
        self.IA = IA
        self.IS = IS
        self.H = H
        self.U = U
        self.R = R         
        return self.ypred

#Auxiliary functions to compute R(t)
#(Fjj - Fii)
    def __prod(self, i, F):
        P = 1
        for j in range(0, len(F)):
            if i != j:
                P = P * (F[j] - F[i])
        return P
##compute g(x)
    def __gx(self, x, F):
        g = 0
        for i in range(len(F)):
            if 0 != self.__prod(i, F): 
                g += np.exp(-F[i]*x)/self.__prod(i, F)
        g = np.prod(F) * g
        return g

#Integral b(t-x)g(x) dx
    def __int(self, b, t, F):
        res = 0
        for x in range(t+1):
            res += b[t - x] * self.__gx(x, F)
        return res

#Compute R(t)
    def Rt(self, cummulativeCases, date, cutoof):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        date = date
        cummulativeCases.index = date 
        cummulativeCases = cummulativeCases[np.argsort(cummulativeCases)]

        #using cummulative cases
        cummulativeCases = np.diff(cummulativeCases[:len(cummulativeCases) + 1])
        #Defining the F matrix array
        try:
            F = np.array([self.kappa, self.gammaA, self.gammaS])
            #initiate a empety list to get result
            res = []
            for t in range(0,len(cummulativeCases)):
                res.append(cummulativeCases[t]/self.__int(cummulativeCases, t, F))
            self.rt = pd.Series(np.array(res))
            self.rt.index  = np.sort(date.iloc[1:])
            idx_start = np.searchsorted(np.cumsum(cummulativeCases),cutoof)
            return(self.rt.iloc[idx_start:])
        except:
            return("Model must be fitted before R(t) could be computed")
        

    def plot(self,local):
        ypred = self.predict(self.x)
        plt.plot(ypred,c='b',label='Predição Infectados')
        plt.plot(self.y,c='r',marker='o', markersize=3,label='Infectados')
        plt.legend(fontsize=15)
        plt.title('Dinâmica do CoviD19 - {}'.format(local),fontsize=20)
        plt.ylabel('Casos COnfirmados',fontsize=15)
        plt.xlabel('Dias',fontsize = 15)
        plt.show()

    def plotDeath(self,local):
        self.predict(self.x)
        plt.plot(self.dpred,c='b',label='Predição mortes')
        plt.plot(self.d,c='r',marker='o', markersize=3,label='mortos')
        plt.legend(fontsize = 15)
        plt.title('Dinâmica do CoviD19 - {}'.format(local),fontsize=20)
        plt.ylabel('Mortos',fontsize=15)
        plt.xlabel('Dias',fontsize=15)
        plt.show()

    def getCoef(self):
        if self.isBetaChange:
            return ['beta1','beta2','dayBetaChange','gammaH','gammaU', 'delta','h','ia0','is0','e0'],[self.beta1,self.beta2,self.dayBetaChange,self.gammaH,self.gammaU,self.delta,self.h,self.ia0,self.is0,self.e0]
        return ['beta','gammaH','gammaU', 'delta','h','ia0','is0','e0'],[self.beta,self.gammaH,self.gammaU,self.delta,self.h,self.ia0,self.is0,self.e0]

    def plotFit(self):
        plt.style.use('seaborn-deep')
        fig, axes = plt.subplots(figsize = (18,8))
        try:
            plt.plot(self.x, self.ypred, label = "Fitted", c = "red")
            plt.scatter(self.x, self.y, label = "Observed", c = "blue")
            plt.legend(loc='upper left')
            plt.show()
        except:
            print("There is no predicted value")
            
    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        #Define empty lists to recive results
        self.__lypred = []
        self.__ldpred = []
        self.__lspred = []
        self.__lepred = []
        self.__lrpred = []
        self.__lhpred = []
        self.__lupred = []
        self.__lIA = []
        self.__lt1 = []
        self.__lIS = []
        self.__lbeta1 = []
        self.__lbeta2 = []
        self.__lgammaH = []
        self.__lgammaU = []
        self.__ldelta = []
        self.__le0 = []
        self.__lia0 = []
        self.__lis0 = []

        self.y = self.y[np.argsort(self.y)]
        self.d = self.d[np.argsort(self.d)]
        casesSeries = self.__genBoot(self.y, times)
        deathSeries = self.__genBoot(self.d, times)
        for i in range(0,len(casesSeries)):
            super().fit(x = range(1,len(self.y) + 1),
                        y = casesSeries[i],
                        d = deathSeries[i])
            super().predict(0)

            self.__lypred.append(self.ypred)
            self.__ldpred.append(self.dpred)
            self.__lhpred.append(self.H)
            self.__lupred.append(self.U)
            self.__lspred.append(self.S)
            self.__lepred.append(self.E)
            self.__lrpred.append(self.R)
            self.__lIA.append(self.IA)
            self.__lIS.append(self.IS)

            self.__lbeta1.append(self.beta1)
            self.__lbeta2.append(self.beta2)
            self.__lgammaH.append(self.gammaH)
            self.__lgammaU.append(self.gammaU)
            self.__ldelta.append(self.delta)
            self.__le0.append(self.e0)
            self.__lia0.append(self.ia0)
            self.__lis0.append(self.is0)
            self.ciResults = {}
        for i, serie in zip([self.__lypred, self.__ldpred, self.__lhpred, self.__lupred, self.__lspred, 
                  self.__lepred, self.__lrpred, self.__lIA, self.__lIS],
                  ["ypred", "dpred", "hpred", "upred", "spred",
                   "epred", "dpred", "IA", "IS"]):
            meanValue, lb, ub = self.__getConfidenceInterval(i,len(i[0]))
            self.ciResults.update({serie: {"meanValue":meanValue,
                 "lb": lb,
                 "ub": ub}})


        #return df with parameters
        self.parameters = pd.DataFrame.from_dict({
            "beta1": self.__lbeta1,
            "beta2": self.__lbeta2,
            "gammaH": self.__lgammaH,
            "gammaU": self.__lgammaU,
            "delta": self.__ldelta,
            "e0": self.__le0,
            "a0": self.__lia0,
            "is0": self.__lis0
        })

    






          
    
      
       
       
      
            


        
        
  
