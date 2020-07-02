#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:18:44 2020
@author: Rafael Veiga rafaelvalenteveiga@gmail.com
@author: matheustorquato matheusft@gmail.com
ADICIONAR OS OUTROS AUTORES
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as spi
import pyswarms as ps
from pyswarms.utils.plotters import plot_cost_history
import pickle as pk
from numbers import Number
import copy
import matplotlib.gridspec as gridspec
from datetime import date, timedelta


class Models:
    def __init__(self,popSize,nCores=None):
        self.isFit=False
        self.isBetaChange = False
        self.isPredict = False
        self.isCI = False
        self.isRT = False
        self.N = popSize
        self.nCores = nCores
    
    def __validadeVar(self,var,name):
        if len(var)<3:
            print('\nthe '+name+' variable has les than 3 elements!\n')
            return False
        for n in var:
            if not isinstance(n, Number):
                print('\nthe elemente '+str(n)+' in '+name+' variable is not numeric!\n')
                return False
        if name=='y':
            flag = 0
            i = 1
            while i < len(var): 
                if(var[i] < var[i - 1]): 
                    flag = 1
                i += 1
            if flag:
                print('\nthe y is not sorted!\n')
                return False
        var = sorted(var) 
        if var[0]<0:
            print('the '+ name + ' can not have negative  '+ str(var[0])+' value!')
            return False
        if name=='y' and var[0]==0:
            print('the y can not have 0 value!')
            return False
        return True
               
        
    def __genBoot(self, series, times = 500):
        series = np.diff(series)
        series = np.insert(series, 0, 1)
        series[series < 0] = 0
        results = []
        for i in range(0,times):
            results.append(np.random.multinomial(n = sum(series), pvals = series/sum(series)))
        return np.array(results)
    
    def __getConfidenceInterval(self, series, level):
        series = np.array(series)
        length = len(series[0])
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
        return (lowerBound, UpperBound)
        
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

class SIR(Models):
    ''' SIR Model'''
    
    def getR0(self):
        if self.isFit:
            if self.isBetaChange:
                return self.beta1/self.gamma
            else:
                return self.beta/self.gamma
        else:
            print("\nThe models is not fitted!\n")
            return None
    
    def __cal_EDO(self,x,beta,gamma):
            t_range = x
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
            t_range = x
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
    
    
    def __residuals(self,coef,x ,y,stand_error):
        if (self.isBetaChange) & (self.dayBetaChange==None):
            S,I,R = self.__cal_EDO_2(x,coef[0],coef[1],coef[2],coef[3])
        elif self.isBetaChange:
            S,I,R = self.__cal_EDO_2(x,coef[0],coef[1],coef[2],self.dayBetaChange)
        else:
            S,I,R = self.__cal_EDO(x,coef[0],coef[1])
        res = (y-(I+R))
        if stand_error:
            res = res / np.sqrt((I+R)+1)
        return res
    
    def __objectiveFunction(self,coef,x ,y,stand_error):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        for i in range(tam2):
            soma[i] = ((self.__residuals(coef[i], x, y, stand_error))**2).mean()
        return soma
    
    def _prepare_bound(self, bound):
        if (self.isBetaChange) & (self.dayBetaChange==None):
            if len(bound[0])==2:
                bound = (bound[0].copy(),bound[1].copy())
                bound[0].append(bound[0][0])
                bound[1].append(bound[1][0])
                bound[0].append(self.x[4])
                bound[1].append(self.x[-5])
                bound[0][3] = self.x[4]
                bound[1][3] = self.x[-5]
            elif len(bound[0]) != 4:
                raise ValueError("Bound of Incorrect size")
        elif self.isBetaChange:
            if len(bound[0])==2:
                bound = (bound[0].copy(),bound[1].copy())
                bound[0].append(bound[0][1])
                bound[1].append(bound[1][1])
            elif len(bound[0]) != 3:
                raise ValueError("Bound of Incorrect size")
        self.bound = bound
    
    def fit(self, y , bound = ([0,1/21],[1,1/5]),stand_error=True, isBetaChange=False, dayBetaChange = None,particles=50,itera=500,c1= 0.5, c2= 0.3, w = 0.9, k=3,norm=1):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
       
        if not self._Models__validadeVar(y,'y'):
            return
        x = range(1,len(y)+1)
        self.isBetaChange = isBetaChange
        self.dayBetaChange = dayBetaChange
        self.y = y
        self.x = x
        self.bound = bound
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm
        
        df = np.array(y)
        self.I0 = df[0]
        self.S0 = 1-self.I0
        self.R0 = 0
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
        optimizer = None
        if bound==None:
            self.bound = bound
            if (isBetaChange) & (dayBetaChange==None):
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=4, options=options)
            elif isBetaChange:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=3, options=options)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=2, options=options)                
        else:
            self._prepare_bound(bound)
            optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=len(bound[0]), options=options,bounds=self.bound)
        cost = pos = None
        if isBetaChange:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=df,stand_error=stand_error,n_processes=self.nCores)
        else:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=df,stand_error=stand_error,n_processes=self.nCores)
        if isBetaChange:
            self.beta1 = pos[0]
            self.gamma = pos[1]
            self.beta2 = pos[2]
            if dayBetaChange==None:
                self.dayBetaChange = pos[3]
            else:
                self.dayBetaChange = dayBetaChange
        else:
            self.beta = pos[0]
            self.gamma = pos[1]
        self.rmse = cost
        self.cost_history = optimizer.cost_history
        self.isFit=True
        
            
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays) 
        self.predictNumDays = numDays
        if self.isBetaChange:
            S,I,R = self.__cal_EDO_2(x,self.beta1,self.gamma,self.beta2,self.dayBetaChange)
        else:
            S,I,R = self.__cal_EDO(x,self.beta,self.gamma)
        self.ypred = I+R
        self.S = S
        self.I = I
        self.R = R
        self.isPredict=True        
        return self.ypred

    def ArangePlots(self,CompartmentPlots):
        
        PlotList=[]
        LabelList=[]
        for i in CompartmentPlots:
        
            if i=='S':
                PlotList.append(self.S)
                LabelList.append('Susceptible individuals')
            # elif i=='E':
            #     PlotList.append(self.E)
            #     LabelList.append('Exposed individuals')
            elif i=='I':
                PlotList.append(self.I)
                LabelList.append('Infected individuals')
            elif i=='R':
                PlotList.append(self.R)
                LabelList.append('Recovered individuals')
            # elif i=='IA':
            #     PlotList.append(self.IA)
            #     LabelList.append('Asymptomatic individuals')
            # elif i=='IS':
            #     PlotList.append(self.IS)
            #     LabelList.append('Symptomatic individuals')
            # elif i=='H':
            #     PlotList.append(self.H)
            #     LabelList.append('Clinic ocupation')
            # elif i=='U':
            #     PlotList.append(self.U)
            #     LabelList.append('ICU ocupation')
            # elif i=='D':
            #     PlotList.append(self.D)
            #     LabelList.append('Cumulative deaths')
            # elif i=='dD':
            #     PlotList.append(np.diff(self.D))
            #     LabelList.append('New deaths')
            elif i=='Y':
                PlotList.append(self.ypred)
                LabelList.append('Cumulative cases')
            elif i=='dY':
                PlotList.append(np.diff(self.ypred))
                LabelList.append('New cases')
            else:
                print('\nThere is no compartment such as "'+str(i)+'" in the model.\n')
               
        return PlotList,LabelList
        

    def plot(self,local=None,InitialDate=None,CompartmentPlots=None,SaveFile=None):
        
        if self.isPredict==False:
            self.predict(0)
        
        
        

        if InitialDate != None:
            initial_date=date(int(InitialDate[0:4]), int(InitialDate[5:7]), int(InitialDate[8:11]))


            dates = []
            dates.append(initial_date.strftime('%Y-%m-%d'))  

            for i in range(len(self.ypred)-1):
                d=initial_date + timedelta(days=i)
                dates.append(d.strftime('%Y-%m-%d'))  

        else:
            dates=np.arange(len(self.ypred))
       
        
        #Plotting
        
        if CompartmentPlots==None:
            
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


            ax.plot(dates,self.ypred,'b-', linewidth=2.5,label='Model')

            ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)


            if self.isBetaChange == True:
            
                ax.axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta change')

         ##################
    
    
    
    
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.65, box.height])
            legend_x = 1
            legend_y = 0.5



        

            ax.tick_params(labelsize=14)
            ax.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

            
                
            ax.set_ylabel('Confirmed cases',fontsize=15)
        
            if InitialDate != None:
                ax.set_xlabel('Days',fontsize=15)
    
    
            ax.xaxis.set_major_locator(plt.MaxNLocator(9))
            plt.setp(ax.get_xticklabels(), rotation=25)


            for tick in ax.get_xticklabels():
                tick.set_fontname("Arial")
            for tick in ax.get_yticklabels():
                tick.set_fontname("Arial")  
        
        
                    
            if local == None:
                fig.suptitle('Model predictions',fontsize=24)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=24)
                
        elif len(CompartmentPlots)==1:
        
            PlotList,LabelList= self.ArangePlots(CompartmentPlots)
            
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


            ax.plot(dates[:len(PlotList[0])],PlotList[0],'b-', linewidth=2.5,label='Model')

            if CompartmentPlots[0]=='Y':
                ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
            elif CompartmentPlots[0]=='dY':
                ax.scatter(dates[:len( self.y)-1], np.diff(self.y),  s=18,color='black',label='Reported data',zorder=3)

            if self.isBetaChange == True:
            
                ax.axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta Change')

         ##################
    
    
    
    
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.65, box.height])
            legend_x = 1
            legend_y = 0.5



        

            ax.tick_params(labelsize=14)
            ax.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

    
                
            ax.set_ylabel(LabelList[0],fontsize=15)
        
            if InitialDate == None:
                ax.set_xlabel('Days',fontsize=15)
    
    
            ax.xaxis.set_major_locator(plt.MaxNLocator(9))
            plt.setp(ax.get_xticklabels(), rotation=25)


            for tick in ax.get_xticklabels():
                tick.set_fontname("Arial")
            for tick in ax.get_yticklabels():
                tick.set_fontname("Arial")  
        
                    
            if local == None:
                fig.suptitle('Model predictions',fontsize=24)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=24)
        
        else:
            
            color=['blue','red','green','darkviolet','orange','darkblue']
            
            
            PlotList,LabelList= self.ArangePlots(CompartmentPlots)
            
            
            if len(CompartmentPlots)==2:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 1, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
            
                ax.append(ax1)
                ax.append(ax2)

            if len(CompartmentPlots)==3:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 2, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,2*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, 1:3]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[1, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)

            if len(CompartmentPlots)==4:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 2, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,2*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, :2])
                ax4 = plt.subplot(gs[1, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
            
            if len(CompartmentPlots)==5:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 3, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,3*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, 1:3]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[1, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, 2:])
                ax4 = plt.subplot(gs[2, :2])
                ax5 = plt.subplot(gs[2, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
                ax.append(ax5)
 

            if len(CompartmentPlots)==6:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 3, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,3*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, :2])
                ax4 = plt.subplot(gs[1, 2:]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax5 = plt.subplot(gs[2, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax6 = plt.subplot(gs[2, 2:])
                
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
                ax.append(ax5)
                ax.append(ax6)
 
  
          
            for k in range(len(ax)):
                
            

                ax[k].grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


                ax[k].plot(dates[:len(PlotList[k])],PlotList[k],color=color[k], linewidth=2.5,label='Model')

                if CompartmentPlots[k]=='Y':
                    ax[k].scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
                elif CompartmentPlots[k]=='dY':
                    ax[k].scatter(dates[:len( self.y)-1], np.diff(self.y),  s=18,color='black',label='Reported data',zorder=3)

                if self.isBetaChange == True:
            
                    ax[k].axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta Change')

             ##################
    


        

                ax[k].tick_params(labelsize=22)
                #ax[k].legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

    
    
                
            
            
                ax[k].set_ylabel(LabelList[k],fontsize=25)
        
                if InitialDate == None:
                    ax[k].set_xlabel('Days',fontsize=25)
    
    
                ax[k].xaxis.set_major_locator(plt.MaxNLocator(9))
                plt.setp(ax[k].get_xticklabels(), rotation=25)


                for tick in ax[k].get_xticklabels():
                    tick.set_fontname("Arial")
                for tick in ax[k].get_yticklabels():
                    tick.set_fontname("Arial")  
        
        
    
            if local == None:
                fig.suptitle('Model predictions',fontsize=35)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=35)
                
            

        if SaveFile != None:
            fig.savefig(SaveFile,bbox_inches='tight')
        
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

    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
            self.lS = self._Models__getConfidenceInterval(self.__bS, level)
            self.lI = self._Models__getConfidenceInterval(self.__bI, level)
            self.lR = self._Models__getConfidenceInterval(self.__bR, level)
        
            if self.isBetaChange:
                self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
                self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
                self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
            else:
                self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
            self.lGamma = self._Models__getConfidenceInterval(self.__bGamma, level)
            
        #Define empty lists to recive results
        self.__bypred = []
        self.__bS = []
        self.__bI = []
        self.__bR = []
        
        if self.isBetaChange:
            self.__bDayBetaChange=[]
            self.__bBeta1 = []
            self.__bBeta2=[]
        else:
            self.__bBeta=[]
            
        self.__bGamma = []
        
        casesSeries = self._Models__genBoot(self.y, times)
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.fit(y = casesSeries[i], bound = self.bound ,stand_error=self.stand_error, isBetaChange=self.isBetaChange,dayBetaChange = self.dayBetaChange,particles=self.particles,itera=self.itera,c1=self.c1,c2= self.c2, w= self.w,k=self.k,norm=self.norm)
            if self.isPredict:
                copia.predict(self.predictNumDays)
            else:
                copia.predict(0)

            self.__bypred.append(copia.ypred)
            self.__bS.append(copia.S)
            self.__bI.append(copia.I)
            self.__bR.append(copia.R)
            if self.isBetaChange:
                self.__bbeta1.append(copia.beta1)
                self.__bbeta2.append(copia.beta2)
                self.__bDayBetaChange.append(copia.dayBetaChange)
            else:
                self.__bbeta.append(copia.beta)
            self.__bgamma.append(copia.gamma)            
        
        self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
        self.lS = self._Models__getConfidenceInterval(self.__bS, level)
        self.lI = self._Models__getConfidenceInterval(self.__bI, level)
        self.lR = self._Models__getConfidenceInterval(self.__bR, level)
        
        if self.isBetaChange:
            self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
            self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
            self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
        else:
            self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
        self.lGamma = self._Models__getConfidenceInterval(self.__bGamma, level)
        self.isCI=True

class SEIR(Models):
    ''' SIR Model'''
    
    def getR0(self):
        if self.isFit:
            if self.isBetaChange:
                return self.beta1/self.gamma
            else:
                return self.beta/self.gamma
        else:
            print("\nThe models is not fitted!\n")
            return None
    
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
    

    def fit(self, y , bound = ([0,1/7,1/6],[1.5,1/4,1/4]) ,stand_error=True, isBetaChange=True,dayBetaChange = None,particles=50,itera=500,c1=0.3,c2= 0.3, w= 0.9,k=3,norm=2):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        
        if not self._Models__validadeVar(y,'y'):
            return
        x = range(1,len(y)+1)
        self.y = y
        dy = np.array(y)
        self.x = x
        self.I0 = np.array(y[0])/self.N
        self.S0 = 1-self.I0
        self.R0 = 0
        self.E0 = 0
        self.mu = 1/(75.51*365)

        self.isBetaChange = isBetaChange
        self.dayBetaChange = dayBetaChange
        
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm


        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
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
        self.bound = bound
        if isBetaChange:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=dy,stand_error=stand_error,n_processes=self.nCores)
        else:
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=dy,stand_error=stand_error,n_processes=self.nCores)
            
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
            
        else:
            self.beta = pos[0]
            self.gamma = pos[1]
            self.sigma = pos[2]
            

        self.rmse = cost
        self.cost_history = optimizer.cost_history
        self.isFit=True
            
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays)
        self.predictNumDays = numDays
        if self.isBetaChange:
            S,E,I,R = self.__cal_EDO_2(x,self.beta1,self.beta2,self.dayBetaChange,self.gamma,self.mu,self.sigma)
        else:
            S,E,I,R = self.__cal_EDO(x,self.beta,self.gamma,self.mu,self.sigma)
        self.ypred = I+R
        self.S = S
        self.E = E
        self.I = I
        self.R = R
        self.isPredict=True        
        return self.ypred

    def ArangePlots(self,CompartmentPlots):
        
        PlotList=[]
        LabelList=[]
        for i in CompartmentPlots:
        
            if i=='S':
                PlotList.append(self.S)
                LabelList.append('Susceptible individuals')
            elif i=='E':
                PlotList.append(self.E)
                LabelList.append('Exposed individuals')
            elif i=='I':
                PlotList.append(self.I)
                LabelList.append('Infected individuals')
            elif i=='R':
                PlotList.append(self.R)
                LabelList.append('Recovered individuals')
            elif i=='Y':
                PlotList.append(self.ypred)
                LabelList.append('Cumulative cases')
            elif i=='dY':
                PlotList.append(np.diff(self.ypred))
                LabelList.append('New cases')
            else:
                print('\nThere is no compartment such as "'+str(i)+'" in the model.\n')
               
        return PlotList,LabelList
        

    def plot(self,local=None,InitialDate=None,CompartmentPlots=None,SaveFile=None):
        
        if self.isPredict==False:
            self.predict(0)
        
        
        

        if InitialDate != None:
            initial_date=date(int(InitialDate[0:4]), int(InitialDate[5:7]), int(InitialDate[8:11]))


            dates = []
            dates.append(initial_date.strftime('%Y-%m-%d'))  

            for i in range(len(self.ypred)-1):
                d=initial_date + timedelta(days=i)
                dates.append(d.strftime('%Y-%m-%d'))  

        else:
            dates=np.arange(len(self.ypred))
       
        
        #Plotting
        
        if CompartmentPlots==None:
            
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


            ax.plot(dates,self.ypred,'b-', linewidth=2.5,label='Model')

            ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)


            if self.isBetaChange == True:
            
                ax.axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta change')

         ##################
    
    
    
    
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.65, box.height])
            legend_x = 1
            legend_y = 0.5



        

            ax.tick_params(labelsize=14)
            ax.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

            
                
            ax.set_ylabel('Confirmed cases',fontsize=15)
        
            if InitialDate != None:
                ax.set_xlabel('Days',fontsize=15)
    
    
            ax.xaxis.set_major_locator(plt.MaxNLocator(9))
            plt.setp(ax.get_xticklabels(), rotation=25)


            for tick in ax.get_xticklabels():
                tick.set_fontname("Arial")
            for tick in ax.get_yticklabels():
                tick.set_fontname("Arial")  
        
        
                    
            if local == None:
                fig.suptitle('Model predictions',fontsize=24)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=24)
                
        elif len(CompartmentPlots)==1:
        
            PlotList,LabelList= self.ArangePlots(CompartmentPlots)
            
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


            ax.plot(dates[:len(PlotList[0])],PlotList[0],'b-', linewidth=2.5,label='Model')

            if CompartmentPlots[0]=='Y':
                ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
            elif CompartmentPlots[0]=='dY':
                ax.scatter(dates[:len( self.y)-1], np.diff(self.y),  s=18,color='black',label='Reported data',zorder=3)

            if self.isBetaChange == True:
            
                ax.axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta Change')

         ##################
    
    
    
    
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.65, box.height])
            legend_x = 1
            legend_y = 0.5



        

            ax.tick_params(labelsize=14)
            ax.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

    
                
            ax.set_ylabel(LabelList[0],fontsize=15)
        
            if InitialDate == None:
                ax.set_xlabel('Days',fontsize=15)
    
    
            ax.xaxis.set_major_locator(plt.MaxNLocator(9))
            plt.setp(ax.get_xticklabels(), rotation=25)


            for tick in ax.get_xticklabels():
                tick.set_fontname("Arial")
            for tick in ax.get_yticklabels():
                tick.set_fontname("Arial")  
        
                    
            if local == None:
                fig.suptitle('Model predictions',fontsize=24)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=24)
        
        else:
            
            color=['blue','red','green','darkviolet','orange','darkblue']
            
            
            PlotList,LabelList= self.ArangePlots(CompartmentPlots)
            
            
            if len(CompartmentPlots)==2:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 1, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
            
                ax.append(ax1)
                ax.append(ax2)

            if len(CompartmentPlots)==3:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 2, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,2*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, 1:3]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[1, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)

            if len(CompartmentPlots)==4:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 2, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,2*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, :2])
                ax4 = plt.subplot(gs[1, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
            
            if len(CompartmentPlots)==5:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 3, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,3*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, 1:3]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[1, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, 2:])
                ax4 = plt.subplot(gs[2, :2])
                ax5 = plt.subplot(gs[2, 2:])
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
                ax.append(ax5)
 

            if len(CompartmentPlots)==6:
           #Criar um grid para as figuras usando GridSpec
                gs = gridspec.GridSpec(nrows = 3, ncols = 4)

        #Definir o tamanho do plot que sera usado para cada plot individula tem o mesmo efieto de quando passado para subplot
                fig=plt.figure(figsize=(2*15.4,3*10))

        #Definir espaco em branco entre os plots
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)

                ax=[]
        #Criar o layout onde os plots serao gerados. E nessa parte que se define o grid
                ax1 = plt.subplot(gs[0, :2]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax2 = plt.subplot(gs[0, 2:])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax3 = plt.subplot(gs[1, :2])
                ax4 = plt.subplot(gs[1, 2:]) #Ininicar um plot em branco no centro da primeira linha (0)
                ax5 = plt.subplot(gs[2, :2])  #Ininicar um plot em branco na primeira posicao da segunda linha
                ax6 = plt.subplot(gs[2, 2:])
                
                
                ax.append(ax1)
                ax.append(ax2)
                ax.append(ax3)
                ax.append(ax4)
                ax.append(ax5)
                ax.append(ax6)
 
  
          
            for k in range(len(ax)):
                
            

                ax[k].grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


                ax[k].plot(dates[:len(PlotList[k])],PlotList[k],color=color[k], linewidth=2.5,label='Model')

                if CompartmentPlots[k]=='Y':
                    ax[k].scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
                elif CompartmentPlots[k]=='dY':
                    ax[k].scatter(dates[:len( self.y)-1], np.diff(self.y),  s=18,color='black',label='Reported data',zorder=3)

                if self.isBetaChange == True:
            
                    ax[k].axvline(self.dayBetaChange, 0, 600,c='r',linestyle='--',label='Beta Change')

             ##################
    


        

                ax[k].tick_params(labelsize=22)
                #ax[k].legend(loc='center left', bbox_to_anchor=(legend_x, legend_y),fontsize=18)

    
    
                
            
            
                ax[k].set_ylabel(LabelList[k],fontsize=25)
        
                if InitialDate == None:
                    ax[k].set_xlabel('Days',fontsize=25)
    
    
                ax[k].xaxis.set_major_locator(plt.MaxNLocator(9))
                plt.setp(ax[k].get_xticklabels(), rotation=25)


                for tick in ax[k].get_xticklabels():
                    tick.set_fontname("Arial")
                for tick in ax[k].get_yticklabels():
                    tick.set_fontname("Arial")  
        
        
    
            if local == None:
                fig.suptitle('Model predictions',fontsize=35)
            else:
                fig.suptitle('Model predictions - '+ local,fontsize=35)
                
            

        if SaveFile != None:
            fig.savefig(SaveFile,bbox_inches='tight')
        
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
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
            self.lS = self._Models__getConfidenceInterval(self.__bS, level)
            self.lE = self._Models__getConfidenceInterval(self.__bE, level)
            self.lI = self._Models__getConfidenceInterval(self.__bI, level)
            self.lR = self._Models__getConfidenceInterval(self.__bR, level)
        
            if self.isBetaChange:
                self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
                self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
                self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
            else:
                self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
            self.lGamma = self._Models__getConfidenceInterval(self.__bGamma, level)
            
        #Define empty lists to recive results
        self.__bypred = []
        self.__bS = []
        self.__bE = []
        self.__bI = []
        self.__bR = []
        
        if self.isBetaChange:
            self.__bDayBetaChange=[]
            self.__bBeta1 = []
            self.__bBeta2=[]
        else:
            self.__bBeta=[]
            
        self.__bGamma = []
        
        casesSeries = self._Models__genBoot(self.y, times)
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.fit(y = casesSeries[i], bound = self.bound ,stand_error=self.stand_error, isBetaChange=self.isBetaChange,dayBetaChange = self.dayBetaChange,particles=self.particles,itera=self.itera,c1=self.c1,c2= self.c2, w= self.w,k=self.k,norm=self.norm)
            if self.isPredict:
                copia.predict(self.predictNumDays)
            else:
                copia.predict(0)

            self.__bypred.append(copia.ypred)
            self.__bS.append(copia.S)
            self.__bE.append(copia.E)
            self.__bI.append(copia.I)
            self.__bR.append(copia.R)
            if self.isBetaChange:
                self.__bbeta1.append(copia.beta1)
                self.__bbeta2.append(copia.beta2)
                self.__bDayBetaChange.append(copia.dayBetaChange)
            else:
                self.__bbeta.append(copia.beta)
            self.__bgamma.append(copia.gamma)            
        
        self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
        self.lS = self._Models__getConfidenceInterval(self.__bS, level)
        self.lE = self._Models__getConfidenceInterval(self.__bE, level)
        self.lI = self._Models__getConfidenceInterval(self.__bI, level)
        self.lR = self._Models__getConfidenceInterval(self.__bR, level)
        
        if self.isBetaChange:
            self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
            self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
            self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
        else:
            self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
        self.lGamma = self._Models__getConfidenceInterval(self.__bGamma, level)
        self.isCI=True
    
class SEIRHUD(Models):
    ''' SEIRHU Model'''
    
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
    
    def __objectiveFunction(self,coef,x ,y,d,hos,u,stand_error):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        if stand_error:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8],coef[i,9])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*self.yWeight +(((d-(D))/np.sqrt(D+1))**2).mean()*self.dWeight
                    soma[i] = (soma[i] + (((hos-(H))/np.sqrt(H+1))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + (((u-(U))/np.sqrt(U+1))**2).mean()*self.uWeight) if u else soma[i]
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange,coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*(1-self.pesoMorte)+(((d-(D))/np.sqrt(D+1))**2).mean()*self.pesoMorte
                    soma[i] = (soma[i] + (((hos-(H))/np.sqrt(H+1))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + (((u-(U))/np.sqrt(U+1))**2).mean()*self.uWeight) if u else soma[i]
            else:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7])
                    soma[i]= (((y-(Nw))/np.sqrt(Nw+1))**2).mean()*(1-self.pesoMorte)+(((d-(D))/np.sqrt(D+1))**2).mean()*self.pesoMorte
                    soma[i] = (soma[i] + (((hos-(H))/np.sqrt(H+1))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + (((u-(U))/np.sqrt(U+1))**2).mean()*self.uWeight) if u else soma[i]
        else:
            if (self.isBetaChange) & (self.dayBetaChange==None):
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8],coef[i,9])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
                    soma[i] = (soma[i] + ((hos-(H))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + ((u-(U))**2).mean()*self.uWeight) if u else soma[i]
            elif self.isBetaChange:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,coef[i,0],coef[i,1],coef[i,2],self.dayBetaChange,coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7],coef[i,8])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
                    soma[i] = (soma[i] + ((hos-(H))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + ((u-(U))**2).mean()*self.uWeight) if u else soma[i]
            else:
                for i in range(tam2):
                    S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,coef[i,0],coef[i,1],coef[i,2],coef[i,3],coef[i,4],coef[i,5],coef[i,6],coef[i,7])
                    soma[i]= ((y-(Nw))**2).mean()*(1-self.pesoMorte)+((d-(D))**2).mean()*self.pesoMorte
                    soma[i] = (soma[i] + ((hos-(H))**2).mean()*self.hosWeight) if hos else soma[i]
                    soma[i] = (soma[i] + ((u-(U))**2).mean()*self.uWeight) if u else soma[i]
        return soma
    def fit(self, y, d, hos=None,u=None,yWeight=1,dWeight = 1,hosWeight=1,uWeight=1, kappa = 1/4,p = 0.2,gammaA = 1/3.5, gammaS = 1/4.001, muH = 0.15,
            muU = 0.4,xi = 0.53,omegaU = 0.29,omegaH=0.14 , bound = [[0,1/8,1/12,0,0],[2,1/4,1/3,0.7,0.35]],
            stand_error = True, isBetaChange = False, dayBetaChange = None, particles = 300, itera = 1000, c1 = 0.1, c2 = 0.3, w = 0.9, k = 5, norm = 2):
        '''
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        
        if not self._Models__validadeVar(y,'y'):
            return
        if not self._Models__validadeVar(d,'d'):
            return
        if hos:
            if not self._Models__validadeVar(hos,'hos'):
                return
        if u:
            if not self._Models__validadeVar(u,'u'):
                return
            
        if len(y)!=len(d):
            print('\ny and d must have the same length\n')
            return
        if hos:
            if len(y)!=len(hos):
                print('\ny and hos must have the same length\n')
                return
        if u:
            if len(y)!=len(u):
                print('\ny and u must have the same length\n')
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
        self.yWeight = yWeight
        self.dWeight = dWeight
        self.hosWeight = hosWeight
        self.uWeight = uWeight
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
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm
        df = np.array(y)
        dd = np.array(d)
        dhos = np.array(hos)
        du = np.array(u)
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
            elif self.isBetaChange:
                if len(bound[0])==8:
                    bound = (bound[0].copy(),bound[1].copy())
                    bound[0].insert(1,bound[0][0])
                    bound[1].insert(1,bound[1][0])
                    
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=9, options=options,bounds=bound)
            else:
                optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=8, options=options,bounds=bound)
                
        cost = pos = None
        self.bound = bound
        #__cal_EDO(self,x,beta,gammaH,gammaU,delta,h,ia0,is0,e0)
        #__cal_EDO_2(self,x,beta1,beta2,tempo,gammaH,gammaU,delta,h,ia0,is0,e0)
        if self.isBetaChange:
            #cost, pos = optimizer.optimize(self.objectiveFunction,itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores, verbose = True)
            cost, pos = optimizer.optimize(self.__objectiveFunction,itera, x = x,y=df,d=dd,hos=dhos,u=du,stand_error=stand_error,n_processes=self.nCores)
        else:
            #cost, pos = optimizer.optimize(self.objectiveFunction, itera, x = x,y=df,d=dd,stand_error=stand_error,n_processes=self.numeroProcessadores, verbose = True)
            cost, pos = optimizer.optimize(self.__objectiveFunction, itera, x = x,y=df,d=dd,hos=dhos,u=du,stand_error=stand_error,n_processes=self.nCores)
            self.beta = pos[0]
            self.gammaH = pos[1]
            self.gammaU = pos[2]
            self.delta = pos[3]
            self.h = pos[4]
            self.ia0 = pos[5]
            self.is0 = pos[6]
            self.e0 = pos[7]
        if self.isBetaChange:
            self.beta1 = pos[0]
            self.beta2 = pos[1]
            
            if self.dayBetaChange==None:
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
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays) 
        self.predictNumDays = numDays
        
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
        self.isPredict=True         
        return self.ypred

#Compute R(t)
    def Rt(self, cutoof):
        #Auxiliary functions to compute R(t)
        #(Fjj - Fii)
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        
        def __prod(i, F):
            P = 1
            for j in range(0, len(F)):
                    if i != j:
                        P = P * (F[j] - F[i])
            return P
        ##compute g(x)
        def __gx(x, F):
            g = 0
            for i in range(len(F)):
                if 0 != self.__prod(i, F): 
                    g += np.exp(-F[i]*x)/__prod(i, F)
            g = np.prod(F) * g
            return g
        #Integral b(t-x)g(x) dx
        def __int( b, t, F):
            res = 0
            for x in range(t+1):
                res += b[t - x] * __gx(x, F)
            return res
        
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        

        cummulativeCases = np.array(self.y)
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
            idx_start = np.searchsorted(np.cumsum(cummulativeCases),cutoof)
            self.isRT=True
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
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
            self.ldpred = self._Models__getConfidenceInterval(self.__bdpred, level)
            self.lH = self._Models__getConfidenceInterval(self.__bH, level)
            self.lS = self._Models__getConfidenceInterval(self.__bS, level)
            self.lE = self._Models__getConfidenceInterval(self.__bE, level)
            self.lR = self._Models__getConfidenceInterval(self.__bR, level)
            self.lH = self._Models__getConfidenceInterval(self.__bH, level)
            self.lU = self._Models__getConfidenceInterval(self.__bU, level)
            self.lIA = self._Models__getConfidenceInterval(self.__bIA, level)
            self.lIS = self._Models__getConfidenceInterval(self.__bIS, level)
        
            if self.isBetaChange:
                self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
                self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
                self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
            else:
                self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
            self.lGammaH = self._Models__getConfidenceInterval(self.__bGammaH, level)
            self.lGammaU = self._Models__getConfidenceInterval(self.__bGammaU, level)
            self.lDelta = self._Models__getConfidenceInterval(self.__bDelta, level)
            self.le0 = self._Models__getConfidenceInterval(self.__be0, level)
            self.lia0 = self._Models.__getConfidenceInterval(self.__bia0, level)
            self.lis0 = self._Models__getConfidenceInterval(self.__bis0, level)
            
        #Define empty lists to recive results
        self.__bypred = []
        self.__bdpred = []
        self.__bS = []
        self.__bE = []
        self.__bR = []
        self.__bH = []
        self.__bU = []
        self.__bIA = []
        self.__bIS = []
        if self.isBetaChange:
            self.__bDayBetaChange=[]
            self.__bBeta1 = []
            self.__bBeta2=[]
        else:
            self.__bBeta=[]
            
        self.__bGammaH = []
        self.__bGammaU = []
        self.__bDelta = []
        self.__be0 = []
        self.__bia0 = []
        self.__bis0 = []

        
        casesSeries = self._Models__genBoot(self.y, times)
        deathSeries = self._Models__genBoot(self.d, times)
        hosSeries = self._Models__genBoot(self.hos,times) if self.hos else None
        uSeries = self._Models.__genBoot(self.u,times) if self.u else None
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.fit(y = casesSeries[i],
                        d = deathSeries[i],
                        hos=hosSeries[i],
                        u = uSeries[i],
                        yWeight=self.yWeight,dWeight = self.dWeight,hosWeight=self.hosWeight,uWeight=self.uWeight,
                        kappa = self.kappa,p = self.p,gammaA = self.gammaA, gammaS = self.gammaS, muH = self.muH,muU = self.muU,xi = self.xi,omegaU = self.omegaU,omegaH=self.omegaH , bound = self.bound,
                        stand_error = self.stand_error, isBetaChange = self.isBetaChange, dayBetaChange = self.dayBetaChange, particles = self.particles, itera = self.itera, c1 = self.c1, c2 = self.c2, w = self.w, k = self.k, norm = self.norm)
            
            if self.isPredict:
                copia.predict(self.predictNumDays)
            else:
                copia.predict(0)
            self.__bypred.append(copia.ypred)
            self.__bdpred.append(copia.dpred)
            self.__bH.append(copia.H)
            self.__bU.append(copia.U)
            self.__bS.append(copia.S)
            self.__bE.append(copia.E)
            self.__bR.append(copia.R)
            self.__bIA.append(copia.IA)
            self.__bIS.append(copia.IS)
            if self.isBetaChange:
                self.__bbeta1.append(copia.beta1)
                self.__bbeta2.append(copia.beta2)
                self.__bDayBetaChange.append(copia.dayBetaChange)
            else:
                self.__bbeta.append(copia.beta)
            self.__bgammaH.append(copia.gammaH)
            self.__bgammaU.append(copia.gammaU)
            self.__bdelta.append(copia.delta)
            self.__be0.append(copia.e0)
            self.__bia0.append(copia.ia0)
            self.__bis0.append(copia.is0)
            
        
        self.lypred = self._Models__getConfidenceInterval(self.__bypred, level)
        self.ldpred = self._Models__getConfidenceInterval(self.__bdpred, level)
        self.lH = self._Models__getConfidenceInterval(self.__bH, level)
        self.lS = self._Models__getConfidenceInterval(self.__bS, level)
        self.lE = self._Models__getConfidenceInterval(self.__bE, level)
        self.lR = self._Models__getConfidenceInterval(self.__bR, level)
        self.lH = self._Models__getConfidenceInterval(self.__bH, level)
        self.lU = self._Models__getConfidenceInterval(self.__bU, level)
        self.lIA = self._Models__getConfidenceInterval(self.__bIA, level)
        self.lIS = self._Models__getConfidenceInterval(self.__bIS, level)
        
        if self.isBetaChange:
            self.lDayBetaChange=self._Models__getConfidenceInterval(self.__bDayBetaChange, level)
            self.lBeta1 = self._Models__getConfidenceInterval(self.__bBeta1, level)
            self.lBeta2 = self._Models__getConfidenceInterval(self.__bBeta2, level)
        else:
            self.lBeta = self._Models__getConfidenceInterval(self.__bBeta, level)
            
        self.lGammaH = self._Models__getConfidenceInterval(self.__bGammaH, level)
        self.lGammaU = self._Models__getConfidenceInterval(self.__bGammaU, level)
        self.lDelta = self._Models__getConfidenceInterval(self.__bDelta, level)
        self.le0 = self._Models__getConfidenceInterval(self.__be0, level)
        self.lia0 = self._Models__getConfidenceInterval(self.__bia0, level)
        self.lis0 = self._Models__getConfidenceInterval(self.__bis0, level)
        
        self.isCI=True