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
import pyswarms.backend as P
from pyswarms.backend.topology import Ring


def load(fileName):
        file = open(fileName,'rb')
        model = pk.load(file)
        return model

class Models:
    def __init__(self,popSize,nCores=None):
        self.isFit=False
        self.BetaChange = 0
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
            aux=0
            for v in var:
                if aux>v:
                    flag=1
                    v1=v
                    v2 =aux
                aux=v
            if flag:
                print('\nthe y is not sorted!\n'+str(v1)+' most be < '+str(v2)+'\n')
                return False
        var = sorted(var) 
        if var[0]<0:
            print('the '+ name + ' can not have negative  '+ str(var[0])+' value!')
            return False
        if name=='y' and var[0]==0:
            print('the y can not have 0 value!')
            return False
        return True
               
    def __changeCases(self, y):
        tam = len(y)
        res = np.ones(tam)
        res[0] = y[0]
        for i in range(1,tam):
            res[i] = y[i]-y[i-1]
        return res
    
    def __genBoot(self, series, times = 500):
        series = np.diff(series)
        series = np.insert(series, 0, 1)
        series[series < 0] = 0
        results = []
        for i in range(0,times):
            results.append(np.random.multinomial(n = sum(series), pvals = series/sum(series)))
        return np.array(results)
    
    def __getConfidenceInterval(self, series, level,isVar):
        series = np.array(series)
        if isVar:
            #Compute mean value
            meanValue = np.mean(series) 
            #Compute deltaStar
            deltaStar = meanValue - series
            #Compute lower and uper bound
            q= (1-level)/2
            deltaL = np.quantile(deltaStar, q = q)
            deltaU = np.quantile(deltaStar, q = 1-q) 

        else:
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
        
    
    
    # def __fitNBeta(self,dim,n_particles,itera,options,objetive_function,BetaChange,bound):
    #     my_topology = Ring()
    #     my_swarm = P.create_swarm(n_particles=n_particles, dimensions=dim, options=options,bounds=bound)
    #     my_swarm.pbest_cost = np.full(n_particles, np.inf)
    #     my_swarm.best_cost = np.inf
        
    #     for i in range(itera):
    #         for a in range(n_particles):
    #             my_swarm.position[a][0:BetaChange] = sorted(my_swarm.position[a][0:BetaChange])
    #             for c in range(1,self.BetaChange):
    #                 if my_swarm.position[a][c-1]+5>=my_swarm.position[a][c]:
    #                     my_swarm.position[a][c]=my_swarm.position[a][c]+5
    #         my_swarm.current_cost = objetive_function(my_swarm.position)
    #         my_swarm.pbest_pos, my_swarm.pbest_cost = P.operators.compute_pbest(my_swarm)
    #         #my_swarm.current_cost[np.isnan(my_swarm.current_cost)]=np.nanmax(my_swarm.current_cost)
    #         #my_swarm.pbest_cost = objetive_function(my_swarm.pbest_pos)
            
            
    #         my_swarm.best_pos, my_swarm.best_cost = my_topology.compute_gbest(my_swarm,options['p'],options['k'])
    #         if i%20==0:
    #             print('Iteration: {} | my_swarm.best_cost: {:.4f} | days: {}'.format(i+1, my_swarm.best_cost, str(my_swarm.pbest_pos[my_swarm.pbest_cost.argmin()])))
    #         my_swarm.velocity = my_topology.compute_velocity(my_swarm,  bounds=bound)
    #         my_swarm.position = my_topology.compute_position(my_swarm,bounds=bound)
    #     final_best_cost = my_swarm.best_cost.copy()
    #     final_best_pos = my_swarm.pbest_pos[my_swarm.pbest_cost.argmin()].copy()
    #     return final_best_pos,final_best_cost
    class Coefficient:
        def __init__(self,Betachange,standardized=True):
            self.dpos = {}
            self.dvar = {}
            self.l_var = []
            self.standardized = standardized
            if standardized:
                self.diff  = []
            self.BetaChange = Betachange            
            
        def addVar(self,var,valor):
            '''
            
            Parameters
            ----------
            var : STRING
                DESCRIPTION.
            valor : number, list(number) or None  
                DESCRIPTION.
    
            Returns
            -------
            None.
    
            '''
            tag = False
            pos = len(self.l_var)
            
            if valor==None:
                self.dpos[var]=pos
                self.dvar[var]=None
                self.l_var.append(valor)
                tag = True
            elif isinstance(valor, Number):
                self.dvar[var] = valor
                tag= True
            elif isinstance(valor, list) or isinstance(valor, tuple):
                if len(valor)==2:
                    if(not isinstance(valor[0], Number) ) or (not isinstance(valor[1], Number)):
                        print('\nthe '+str(var)+' variable value is not a number: '+str(valor)+' !\n')
                    elif valor[0]<valor[1]:
                        self.dpos[var] = pos
                        self.dvar[var] = None
                        self.l_var.append(valor)
                        if self.standardized:
                            self.diff.append(valor[1]-valor[0])
                        tag = True
                    else:
                        print('\nthe '+str(var)+' variable value is in wrong order: '+str(valor)+' !\n') 
            else:
                 print('\nthe '+str(var)+' variable value is in wrong format: '+str(valor)+' !\n')   
            return tag
            
        def  getDimention(self):
            return len(self.l_var)
        
        def getDic(self,coef):
            out = {}
            if self.standardized==False:
                for var in self.dvar.keys():
                    if self.dvar[var]!=None:
                        out[var] = self.dvar[var]
                    else:
                        pos = self.dpos[var]
                        out[var]=coef[pos]
            else:
                for var in self.dvar.keys():
                    if self.dvar[var]!=None:
                        out[var] = self.dvar[var]
                    else:
                        pos = self.dpos[var]
                        out[var] = self.diff[pos]*coef[pos] + self.l_var[pos][0]
            beta = []
            tam = self.BetaChange+1
            for i in range(tam):
                nome = 'beta_'+str(i)
                beta.append(out[nome])
                del out[nome]
            out['beta']=beta
            day = []
            tam=tam-1
            for i in range(tam):
                nome = 'dayBetaChange_'+str(i)
                day.append(out[nome])
                del out[nome]
            out['dayBetaChange']=day
                
            return out
        
        def addBeta(self,beta):
            tag=True
            if beta==None:
                beta = []
                for i in range(self.BetaChange+1):
                    beta.append([0,2])
            if isinstance(beta, Number):
                beta = [beta]
            if (len(beta)-1)==self.BetaChange:
               for i in range(self.BetaChange+1):
                   if self.addVar('beta_'+str(i),beta[i])==False:
                       tag=False
            else: 
               print('beta with wrong value: ' + str(beta) + ' !\n')
               tag=False
            return tag
   
        def addDay(self,dayBetaChange):
            if dayBetaChange==None:
                dayBetaChange = [dayBetaChange]
            if isinstance(dayBetaChange, Number):
                dayBetaChange = [dayBetaChange]
                
            tag=True
            if len(dayBetaChange)==self.BetaChange:
               for i in range(self.BetaChange):
                   if self.addVar('dayBetaChange_'+str(i),dayBetaChange[i])==False:
                       tag=False
            else:
                if dayBetaChange!=[None]:    
                    print('dayBetaChange with wrong value: ' + str(dayBetaChange) + ' !\n')
                    tag = False
            return tag
        
        def getBound(self):
            dim = self.getDimention()
            if dim==0:
                bound=None
            else:
                bound = ([],[])
            if self.standardized:
                for i in range(dim):
                    bound[0].append(0)
                    bound[1].append(1)
            else:
                for i in range(dim):
                    bound[0].append(self.l_var[i][0])
                    bound[1].append(self.l_var[i][1])
                               
            return bound
        
            

class SIR(Models):
    ''' SIR Model'''
    
    def getR0(self):
        if self.isFit:
            if self.BetaChange>0:
                return self.beta[0]/self.gamma
            else:
                return self.beta/self.gamma
        else:
            print("\nThe models is not fitted!\n")
            return None
    
    def __cal_EDO(self,x,beta,gamma,I0,R0):
            t_range = x
            def SIR_diff_eqs(INP, t, beta, gamma):
                Y = np.zeros((3))
                V = INP
                Y[0] = - beta * V[0] * V[1]                 #S
                Y[1] = beta * V[0] * V[1] - gamma * V[1]    #I
                Y[2] = gamma * V[1]                         #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (1-(I0+R0), I0,R0), t_range,args=(beta, gamma))
            
            S=result_fit[:, 0]*self.N
            R=result_fit[:, 2]*self.N
            I=result_fit[:, 1]*self.N
            
            return S,I,R
        
    def __cal_EDO_2(self,x,tVar,betaVar,gamma,I0,R0):
            t_range = x
            def beta(t,tVar,betaVar):
                for i in range(len(tVar)):
                    if t<tVar[i]:
                        return betaVar[i]
                return betaVar[i+1]

            
            def SIR_diff_eqs(INP, t, tVar,betaVar, gamma):
                Y = np.zeros((3))
                V = INP
                Y[0] = - beta(t,tVar,betaVar) * V[0] * V[1]                 #S
                Y[1] = beta(t,tVar,betaVar) * V[0] * V[1] - gamma * V[1]    #I
                Y[2] = gamma * V[1]                         #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, ( 1-(I0+R0), I0,R0), t_range,
                                    args=(tVar, betaVar,gamma))
            
            S=result_fit[:, 0]*self.N
            R=result_fit[:, 2]*self.N
            I=result_fit[:, 1]*self.N
            
            return S,I,R
    
    
    def _residuals(self,coef):
        dic = self.coef.getDic(coef)
        dayBetaChange = dic['dayBetaChange']
        beta = dic['beta']
        gamma = dic['gamma']
        I0 = dic['I0']
        R0 = dic['I0']
        if self.BetaChange>0:        
            S,I,R = self.__cal_EDO_2(self.x,dayBetaChange,beta,gamma,I0,R0)
        else:
            S,I,R = self.__cal_EDO(self.x,beta[0],gamma,I0,R)
        aux = I+R
        if self.fittingByCumulativeCases:
            aux2 = self.y-aux 
        else:
            aux = self._Models__changeCases(aux) 
            aux2 = self.NC - aux
        
        if self.stand_error:
            aux2 = aux2 / np.sqrt(aux+1)
        return aux2
    
    def _objectiveFunction(self,coef):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        for i in range(tam2):
            soma[i] = ((self._residuals(coef[i]))**2).mean()
        return soma
    
    
    def __fit(self,x, y ,c1= 0.6, c2= 0.5, w = 0.9):
        options = {'c1': c1, 'c2': c2, 'w': w}
        self.x = x
        self.y = y
        self.NC = self._Models__changeCases(self.y)
        
        dim = self.coef.getDimention()
        bound = self.coef.getBound()
        par=ps.backend.generators.create_swarm(10,dimensions=dim,bounds=bound)
        for i in range(5):
            par.position[i]=self.pos
        optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=dim, options=options,bounds=bound,init_pos=par.position)
        cost, pos = optimizer.optimize(self._objectiveFunction, 1500,n_processes=self.nCores)
        dic = self.coef.getDic(pos)
        self.beta = dic['beta']
        if self.BetaChange>0:
            self.dayBetaChange = dic['dayBetaChange']
        self.gamma = dic['gamma']
        self.S0=dic['S0']
        self.I0=dic['I0']
        self.R0=dic['R0']
        self.predict(0)    

    
    def fit(self,x, y ,fittingByCumulativeCases=True,beta = None,dayBetaChange=None,gamma=[1/21,1/5],I0=None,R0=None,stand_error=True, BetaChange=0,particles=100,itera=1000,c1= 0.5, c2= 0.3, w = 0.9, k=3,norm=2,standarPSO=True):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        if I0==None:
            I0=self.y[0] / self.N
        if R0==None:
           R0=0      
        self.coef = self.Coefficient(BetaChange,standarPSO)
        self.BetaChange = BetaChange
        if not self.coef.addBeta(beta):
                return
        if self.BetaChange!=0:
            if dayBetaChange==None:
                dayBetaChange = []
                for i in range(BetaChange):
                    dayBetaChange.append([])
                    dayBetaChange[i].append(x[5])
                    dayBetaChange[i].append(x[-4])
                if not self.coef.addDay(dayBetaChange):
                    return
            else:
                if not self.coef.addDay(dayBetaChange):
                    return
            if not self.coef.addDay(dayBetaChange):
                return
        if not self.coef.addVar('gamma',gamma):
            return
        if not self.coef.addVar('I0',I0):
            return
        if not self.coef.addVar('R0',R0):
            return
        if not self._Models__validadeVar(y,'y'):
            return
        self.y = np.array(y)
        self.x = np.array(x)
        self.fittingByCumulativeCases = fittingByCumulativeCases
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm
        self.standarPSO = standarPSO
        self.NC = self._Models__changeCases(self.y)
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
        dim = self.coef.getDimention()
        
        if dim==0:
            self.beta = beta
            self.dayBetaChange = dayBetaChange
            self.gamma = gamma
            self.I0=I0
            self.R0=R0
        else:
            optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=dim, options=options,bounds=self.coef.getBound())
            cost, pos = optimizer.optimize(self._objectiveFunction, itera,n_processes=self.nCores)
            self.pos = pos
            dic = self.coef.getDic(pos)
            self.beta = dic['beta']
            if self.BetaChange>0:
                self.dayBetaChange = dic['dayBetaChange']
            self.gamma = dic['gamma']
            self.I0=dic['I0']
            self.R0=dic['R0']
            self.cost_history = optimizer.cost_history
        self.isFit=True
        self.predict(0)
       
            
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        
        x = np.arange(self.x[0], self.x[-1] + 1+numDays) 
        self.predictNumDays = numDays
        if self.BetaChange==0:
            S,I,R = self.__cal_EDO(x,self.beta[0],self.gamma,self.I0,self.R0)
        else:
            S,I,R = self.__cal_EDO_2(x,self.dayBetaChange,self.beta,self.gamma,self.I0,self.R0)
        self.ypred = I+R
        self.xpred = x
        self.S = S
        self.I = I
        self.R = R
        self.NCpred =self._Models__changeCases(self.ypred)
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
            elif i=='NC':
                PlotList.append(self.NCpred)
                LabelList.append('New cases')
            else:
                print('\nThere is no compartment such as "'+str(i)+'" in the model.\n')            
        return PlotList,LabelList
        

    def plot(self,local=None,InitialDate=None,CompartmentPlots=None,SaveFile=None):      

        if InitialDate != None:
            initial_date=date(int(InitialDate[0:4]), int(InitialDate[5:7]), int(InitialDate[8:11]))
            dates = []
            dates.append(initial_date.strftime('%Y-%m-%d'))  
            for i in range(len(self.ypred)-1):
                d=initial_date + timedelta(days=i)
                dates.append(d.strftime('%Y-%m-%d'))  
        else:
            dates=self.xpred
        #Plotting
        
        if CompartmentPlots==None:
            
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)


            ax.plot(dates,self.ypred,'b-', linewidth=2.5,label='Model')

            ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
            if self.BetaChange >0:
                for day in self.dayBetaChange:
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')

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

            if self.BetaChange >0:
                for day in self.dayBetaChange: 
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta Change')

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

                if self.BetaChange > 0:
                    for day in self.dayBetaChange:
                        ax[k].axvline(day, 0, 600,c='r',linestyle='--',label='Beta Change')

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

    
    def getEstimation(self):
        return dict(zip(['S','I','R','Cumulative_cases_predict','new_cases_predict'],[self.S,self.I,self.R,self.ypred,self.NCpred]))
    
    def getCoef(self):
        res = {}
        if self.isFit:
            if self.BetaChange==0:
                res['beta']=self.beta[0]
            else:
                for i in range(len(self.beta)):
                    res['beta'+str(i)] = self.beta[i]
                for i in range(len(self.dayBetaChange)):
                    res['dayBetaChange'+str(i)] = self.dayBetaChange[i]
            res['gamma'] = self.gamma
            res['I0'] = self.I0
            res['R0'] = self.R0
        else:
            print('The model is not fitted!\n')
            res = None
        return res
                    

    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        
        if self.coef.getDimention()==0:
            print('none of coefficients are fitted\n')
            return
        
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
            self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
            self.lI = self._Models__getConfidenceInterval(self._bI, level,False)
            self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
            self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
            self.lCoef = {}
            for k in self.bcoef.keys():
                self.lCoef[k]=self._Models__getConfidenceInterval(self._bCoef[k], level,True)
                        
        #Define empty lists to recive results
        self._bypred = []
        self._bS = []
        self._bI = []
        self._bR = []
        self._bNCpred = []
        listCoef = tuple(self.getCoef().keys())
        self._bCoef={}
        for k in listCoef:
            self._bCoef[k]=[]
        
        casesSeries = self._Models__genBoot(self.y, times)
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.__fit(x=self.x, y=casesSeries[i])
            copiaCoef = copia.getCoef()
            for k in listCoef:
                self._bCoef[k].append(copiaCoef[k])
            self._bypred.append(copia.ypred)
            self._bS.append(copia.S)
            self._bI.append(copia.I)
            self._bR.append(copia.R)
            self._bNCpred.append(copia.NCpred)
                       
        
        self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
        self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
        self.lI = self._Models__getConfidenceInterval(self._bI, level,False)
        self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
        self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
        self.lCoef = {}
        
        for k in listCoef:
            self.lCoef[k] = self._Models__getConfidenceInterval(self._bCoef[k], level,True)
        self.isCI=True

class SEIR(Models):
    ''' SIR Model'''
    
    def getR0(self):
        if self.isFit:
            if self.BetaChange>0:
                return self.beta[0]/self.gamma
            else:
                return self.beta/self.gamma
        else:
            print("\nThe models is not fitted!\n")
            return None
    
    def __cal_EDO(self,x,beta,gamma,mu,sigma,E0,I0,R0):
            t_range = x
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
            result_fit = spi.odeint(SEIR_diff_eqs, (1-(E0+I0+R0),E0,I0,R0), t_range,
                                    args=(beta, gamma,mu,sigma))
            
            S=result_fit[:, 0]*self.N
            E=result_fit[:, 1]*self.N
            I=result_fit[:, 2]*self.N
            R=result_fit[:, 3]*self.N
            
            return S,E,I,R
      
    def __cal_EDO_2(self,x,tVar,betaVar,gamma,mu,sigma,E0,I0,R0):
            t_range = x
            def beta(t,tVar,betaVar):
                for i in range(len(tVar)):
                    if t<tVar[i]:
                        return betaVar[i]
                return betaVar[i+1]
            
            def SEIR_diff_eqs(INP, t, tVar,betaVar, gamma,mu,sigma):
                Y = np.zeros((4))
                V = INP
                Y[0] = mu - beta(t,tVar,betaVar) * V[0] * V[2] - mu * V[0]  # Susceptile
                Y[1] = beta(t,tVar,betaVar) * V[0] * V[2] - sigma * V[1] - mu * V[1] # Exposed
                Y[2] = sigma * V[1] - gamma * V[2] - mu * V[2] # Infectious
                Y[3] = gamma * V[2] #recuperado
                return Y   # For odeint

                return Y
            result_fit = spi.odeint(SEIR_diff_eqs, (1-(E0+I0+R0),E0,I0,R0), t_range,
                                    args=(tVar, betaVar, gamma,mu,sigma))
            
            S=result_fit[:, 0]*self.N
            E=result_fit[:, 1]*self.N
            I=result_fit[:, 2]*self.N
            R=result_fit[:, 3]*self.N
            
            return S,E,I,R
        
    def _residuals(self,coef):
        dic = self.coef.getDic(coef)
        dayBetaChange = dic['dayBetaChange']
        beta = dic['beta']
        gamma = dic['gamma']
        mu = dic['mu']
        sigma = dic['sigma']
        E0 = dic['E0']
        I0 = dic['I0']
        R0 = dic['R0']
        if self.BetaChange>0:        
            S,E,I,R = self.__cal_EDO_2(self.x,dayBetaChange,beta,gamma,mu,sigma,E0,I0,R0)
        else:
            S,E,I,R= self.__cal_EDO(self.x,beta[0],gamma,mu,sigma,E0,I0,R0)
        aux = I+R
        if self.fittingByCumulativeCases:
            aux2 = self.y-aux 
        else:
            aux = self._Models__changeCases(aux) 
            aux2 = self.NC - aux
        
        if self.stand_error:
            aux2 = aux2 / np.sqrt(aux+1)
        return aux2    
    
    def _objectiveFunction(self,coef):
        tam2 = len(coef[:,0])
        soma = np.zeros(tam2)
        for i in range(tam2):
            soma[i] = ((self._residuals(coef[i]))**2).mean()
        return soma
    
    def __fit(self,x, y ,c1= 0.6, c2= 0.5, w = 0.9):
        options = {'c1': c1, 'c2': c2, 'w': w}
        self.x = x
        self.y = y
        self.NC = self._Models__changeCases(self.y)
        dim = self.coef.getDimention()
        bound = self.coef.getBound()
        par=ps.backend.generators.create_swarm(10,dimensions=dim,bounds=bound)
        for i in range(5):
            par.position[i]=self.pos
        optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=dim, options=options,bounds=bound,init_pos=par.position)
        cost, pos = optimizer.optimize(self._objectiveFunction, 1500,n_processes=self.nCores)
        dic = self.coef.getDic(pos)
        self.beta = dic['beta']
        if self.BetaChange>0:
            self.dayBetaChange = dic['dayBetaChange']
        self.gamma = dic['gamma']
        self.mu = dic['mu']
        self.sigma = dic['sigma']
        self.E0 = dic['E0']
        self.I0 = dic['I0']
        self.R0 = dic['R0']
        self.predict(0)    
    
    def fit(self,x, y ,fittingByCumulativeCases=True,beta=None,dayBetaChange=None,gamma=[1/21,1/5],mu=1/(75.51*365),sigma=[1/6,1/4],E0=None,I0=None,R0=None,stand_error=True,BetaChange=0 ,particles=100,itera=1000,c1=0.5,c2= 0.3, w= 0.9,k=3,norm=2,standarPSO=True):
        '''
        x = dias passados do dia inicial 1
        y = numero de casos
        bound = intervalo de limite para procura de cada parametro, onde None = sem limite
        
        bound => (lista_min_bound, lista_max_bound)
        '''
        if I0==None:
            I0=self.y[0] / self.N
        if E0==None:
           E0=0
        if R0==None:
           R0=0
        self.coef = self.Coefficient(BetaChange,standarPSO)
        self.BetaChange = BetaChange
        if not self.coef.addBeta(beta):
                return
        if self.BetaChange!=0:
            if dayBetaChange==None:
                dayBetaChange = []
                for i in range(BetaChange):
                    dayBetaChange.append([])
                    dayBetaChange[i].append(x[5])
                    dayBetaChange[i].append(x[-4])
                if not self.coef.addDay(dayBetaChange):
                    return
            else:
                if not self.coef.addDay(dayBetaChange):
                    return
            if not self.coef.addDay(dayBetaChange):
                return
        if not self.coef.addVar('gamma',gamma):
            return
        if not self.coef.addVar('mu',mu):
            return
        if not self.coef.addVar('sigma',sigma):
            return
        if not self.coef.addVar('E0',E0):
            return
        if not self.coef.addVar('I0',I0):
            return
        if not self.coef.addVar('R0',R0):
            return
        if not self._Models__validadeVar(y,'y'):
            return

        self.y = np.array(y)
        self.x = np.array(x)
        self.fittingByCumulativeCases = fittingByCumulativeCases
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm
        self.standarPSO = standarPSO
        self.NC = self._Models__changeCases(self.y)
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
        dim = self.coef.getDimention()

        if dim==0:
            self.beta = beta
            self.dayBetaChange = dayBetaChange
            self.gamma = gamma
            self.mu = mu
            self.sigma = sigma
            self.E0 = E0
            self.I0 = I0
            self.R0 = R0
        else:
            optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=dim, options=options,bounds=self.coef.getBound())
            cost, pos = optimizer.optimize(self._objectiveFunction, itera,n_processes=self.nCores)
            self.pos = pos
            dic = self.coef.getDic(pos)
            self.beta = dic['beta']
            if self.BetaChange>0:
                self.dayBetaChange = dic['dayBetaChange']
            self.gamma = dic['gamma']
            self.mu = dic['mu']
            self.sigma = dic['sigma']
            self.E0 = dic['E0']
            self.I0 = dic['I0']
            self.R0 = dic['R0']
            self.cost_history = optimizer.cost_history
        self.isFit=True
        self.predict(0)
        
    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        
        x = np.arange(self.x[0], self.x[-1] + 1+numDays) 
        self.predictNumDays = numDays
        if self.BetaChange==0:
            S,E,I,R = self.__cal_EDO(x,self.beta[0],self.gamma,self.mu,self.sigma,self.E0,self.I0,self.R0)
        else:
            S,E,I,R = self.__cal_EDO_2(x,self.dayBetaChange,self.beta,self.gamma,self.mu,self.sigma,self.E0,self.I0,self.R0)
        self.ypred = I+R
        self.xpred = x
        self.S = S
        self.E = E
        self.I = I
        self.R = R
        self.NCpred =self._Models__changeCases(self.ypred)
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
            elif i=='NC':
                PlotList.append(self.NCpred)
                LabelList.append('New cases')
            else:
                print('\nThere is no compartment such as "'+str(i)+'" in the model.\n')
               
        return PlotList,LabelList
        

    def plot(self,local=None,InitialDate=None,CompartmentPlots=None,SaveFile=None):
        
        if InitialDate != None:
            initial_date=date(int(InitialDate[0:4]), int(InitialDate[5:7]), int(InitialDate[8:11]))
            dates = []
            dates.append(initial_date.strftime('%Y-%m-%d'))  
            for i in range(len(self.ypred)-1):
                d=initial_date + timedelta(days=i)
                dates.append(d.strftime('%Y-%m-%d')) 
        else:
            dates=self.xpred
       
        
        #Plotting
        
        if CompartmentPlots==None:           
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)
            ax.plot(dates,self.ypred,'b-', linewidth=2.5,label='Model')
            ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
            if self.BetaChange >0:
                for day in self.dayBetaChange:
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')
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
            if self.BetaChange >0:
                for day in self.dayBetaChange:
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')
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
                gs = gridspec.GridSpec(nrows = 1, ncols = 4)
                fig=plt.figure(figsize=(2*15.4,10))
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)
                ax=[]
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

                if self.BetaChange >0:
                    for day in self.dayBetaChange:
                        ax[k].axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')

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
        res = {}
        if self.isFit:
            if self.BetaChange==0:
                res['beta']=self.beta[0]
            else:
                for i in range(len(self.beta)):
                    res['beta'+str(i)] = self.beta[i]
                for i in range(len(self.dayBetaChange)):
                    res['dayBetaChange'+str(i)] = self.dayBetaChange[i]
            res['gamma'] = self.gamma
            res['mu'] = self.mu
            res['sigma'] = self.sigma
            res['E0'] = self.E0
            res['I0'] = self.I0
            res['R0'] = self.R0
        else:
            print('The model is not fitted!\n')
            res = None
        return res
    
    
    def getEstimation(self):
        return dict(zip(['S','E','I','R','Cumulative_cases_predict','new_cases_predict'],[self.S,self.E,self.I,self.R,self.ypred,self.NCpred]))    
    
    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        if self.coef.getDimention()==0:
            print('none of coefficients are fitted\n')
            return
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
            self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
            self.lE = self._Models__getConfidenceInterval(self._bE, level,False)
            self.lI = self._Models__getConfidenceInterval(self._bI, level,False)
            self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
            self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
            self.lCoef = {}
            for k in self.bcoef.keys():
                self.lCoef[k]=self._Models__getConfidenceInterval(self._bCoef[k], level,True)
                        
        #Define empty lists to recive results
        self._bypred = []
        self._bS = []
        self._bE = []
        self._bI = []
        self._bR = []
        self._bNCpred = []
        listCoef = tuple(self.getCoef().keys())
        self._bCoef={}
        for k in listCoef:
            self._bCoef[k]=[]
        casesSeries = self._Models__genBoot(self.y, times)
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.__fit(x=self.x, y=casesSeries[i])
            copiaCoef = copia.getCoef()
            for k in listCoef:
                self._bCoef[k].append(copiaCoef[k])
            self._bypred.append(copia.ypred)
            self._bS.append(copia.S)
            self._bE.append(copia.E)
            self._bI.append(copia.I)
            self._bR.append(copia.R)
            self._bNCpred.append(copia.NCpred)
                       
        
        self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
        self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
        self.lE = self._Models__getConfidenceInterval(self._bE, level,False)
        self.lI = self._Models__getConfidenceInterval(self._bI, level,False)
        self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
        self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
        self.lCoef = {}
        
        for k in listCoef:
            self.lCoef[k] = self._Models__getConfidenceInterval(self._bCoef[k], level,True)
        self.isCI=True
    
class SEIIHURD(Models):
    ''' SEIRHU Model'''
    
    def __cal_EDO(self,x,beta,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH,E0,Ia0,Is0,H0,U0,R0,D0):
            t_range = x
            def SIR_diff_eqs(INP, t, beta,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH):
                Y = np.zeros((9))
                V = INP
                Y[0] = - beta*V[0]*(V[3] + delta*V[2])                    #S
                Y[1] = beta*V[0]*(V[3] + delta*V[2]) -kappa * V[1]
                Y[2] = (1-p)*kappa*V[1] - gammaA*V[2]
                Y[3] = p*kappa*V[1] - gammaS*V[3]
                Y[4] = h*xi*gammaS*V[3] + (1-muU + omegaU*muU)*gammaU*V[5] -gammaH*V[4]
                Y[5] = h*(1-xi)*gammaS*V[3] +omegaH*gammaH*V[4] -gammaU*V[5]
                Y[6] = gammaA*V[2] + (1-(muH))*(1-omegaH)*gammaH*V[4] + (1-h)*gammaS*V[3]
                Y[7] = (1-omegaH)*muH*gammaH*V[4] + (1-omegaU)*muU*gammaU*V[5]#R
                Y[8] = p*kappa*V[1] 
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (1-(E0+Ia0+Is0),E0,Ia0,Is0,H0,U0,R0,D0,Ia0+Is0), t_range,
                                    args=(beta,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH))
            
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
        
    def __cal_EDO_2(self,x,tVar,betaVar,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH,E0,Ia0,Is0,H0,U0,R0,D0):
            t_range = x
            def beta(t,tVar,betaVar):
                for i in range(len(tVar)):
                    if t<tVar[i]:
                        return betaVar[i]
                return betaVar[i+1]

            delta = np.array(delta)
            def SIR_diff_eqs(INP, t, tVar,betaVar,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH):
                
                Y = np.zeros((9))
                V = INP
                Y[0] = - beta(t,tVar,betaVar)*V[0]*(V[3] + delta*V[2])                    #S
                Y[1] = beta(t,tVar,betaVar)*V[0]*(V[3] + delta*V[2]) -kappa * V[1]
                Y[2] = (1-p)*kappa*V[1] - gammaA*V[2]
                Y[3] = p*kappa*V[1] - gammaS*V[3]
                Y[4] = h*xi*gammaS*V[3] + (1-muU + omegaU*muU)*gammaU*V[5] -gammaH*V[4]
                Y[5] = h*(1-xi)*gammaS*V[3] +omegaH*gammaH*V[4] -gammaU*V[5]
                Y[6] = gammaA*V[2] + (1-(muH))*(1-omegaH)*gammaH*V[4] + (1-h)*gammaS*V[3]
                Y[7] = (1-omegaH)*muH*gammaH*V[4] + (1-omegaU)*muU*gammaU*V[5]#R
                Y[8] = p*kappa*V[1]                      #R
                
                return Y
            result_fit = spi.odeint(SIR_diff_eqs, (1-(E0+Ia0+Is0),E0,Ia0,Is0,H0,U0,R0,D0,Ia0+Is0), t_range,
                                    args=(tVar,betaVar,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH))
            
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
    
    def _objectiveFunction(self,coef):
        tam2 = len(coef[:,0])
        res = np.zeros(tam2)
        for i in range(tam2):
            res[i] = self._residuals(coef[i])
        return res
    
    def _residuals(self,coef):
        dic = self.coef.getDic(coef)
        dayBetaChange = dic['dayBetaChange']
        beta = dic['beta']
        gammaH = dic['gammaH']
        gammaU = dic['gammaU']
        delta = dic['delta']
        kappa = dic['kappa']
        h = dic['h']
        p = dic['p']
        gammaA = dic['gammaA']
        gammaS = dic['gammaS']
        muH = dic['muH']
        muU = dic['muU']
        xi = dic['xi']
        omegaU = dic['omegaU']
        omegaH = dic['omegaH']
        E0 = dic['E0']
        Ia0 = dic['Ia0']
        Is0 = dic['Is0']
        H0 = dic['H0']
        U0 = dic['U0']
        R0 = dic['R0']
        D0 = dic['D0']
        
        if self.BetaChange>0:        
            S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(self.x,dayBetaChange,beta,gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH,E0,Ia0,Is0,H0,U0,R0,D0)
        else:
            S,E,IA,IS,H,U,R,D,Nw= self.__cal_EDO(self.x,beta[0],gammaH,gammaU,delta,kappa,h,p,gammaA,gammaS,muH,muU,xi,omegaU,omegaH,E0,Ia0,Is0,H0,U0,R0,D0)
        
        if self.fittingByCumulativeCases:
            auxy = self.y-Nw 
        else:
            aux = self._Models__changeCases(Nw) 
            auxy = self.NC - aux
        auxd = self.d-D
        auxhos = 0
        if self.hos:
            auxhos=self.hos - H
            if self.stand_error:
                auxhos =((auxhos / np.sqrt(H+1))**2).mean()*self.hosWeight
            else:
                auxhos =(auxhos**2).mean()*self.hosWeight
        auxu = 0
        if self.u:
            auxu=self.u-U
            if self.stand_error:
                auxu = ((auxu / np.sqrt(U+1))**2).mean()*self.uWeight
            else:
                auxu = (auxu**2).mean()*self.uWeight
                
        if self.stand_error:
            auxy = ((auxy / np.sqrt(Nw+1))**2).mean()*self.yWeight
            auxd = ((auxd / np.sqrt(D+1))**2).mean()*self.dWeight
            
        else:
            auxy = (auxy**2).mean()*self.yWeight
            auxd = (auxd**2).mean()*self.dWeight
            
            
        return auxy + auxd + auxhos + auxu
    
    def __fit(self,x, y ,d,hos=None,u=None,c1= 0.6, c2= 0.5, w = 0.9):
        options = {'c1': c1, 'c2': c2, 'w': w}
        self.x = x
        self.y = y
        self.d = d
        self.hos=hos
        self.u = u
        self.NC = self._Models__changeCases(self.y)
        
        dim = self.coef.getDimention()
        bound = self.coef.getBound()
        par=ps.backend.generators.create_swarm(10,dimensions=dim,bounds=bound)
        for i in range(5):
            par.position[i]=self.pos
        optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=dim, options=options,bounds=bound,init_pos=par.position)
        cost, pos = optimizer.optimize(self._objectiveFunction, 1500,n_processes=self.nCores)
        dic = self.coef.getDic(pos)
        self.beta = dic['beta']
        if self.BetaChange>0:
            self.dayBetaChange = dic['dayBetaChange']
        self.gammaH = dic['gammaH']
        self.gammaU = dic['gammaU']
        self.delta = dic['delta']
        self.kappa = dic['kappa']
        self.h = dic['h']
        self.p=dic['p']
        self.gammaA = dic['gammaA']
        self.gammaS = dic['gammaS']
        self.muH = dic['muH']
        self.muU = dic['muU']
        self.xi = dic['xi']
        self.omegaU = dic['omegaU']
        self.omegaH = dic['omegaH']
        self.E0 = dic['E0']
        self.Ia0 = dic['Ia0']
        self.Is0 = dic['Is0']
        self.H0 = dic['H0']
        self.U0 = dic['U0']
        self.R0 = dic['R0']
        self.D0 = dic['D0']
        self.predict(0)
    
    def fit(self, x, y, d,hos=None,u=None,beta=None,dayBetaChange=None,BetaChange = 0,gammaH=[1/8,1/4],gammaU=[1/12,1/3],delta=[0,0.7],kappa = 1/4,h=[0,0.35],p = 0.2,gammaA = 1/3.5,gammaS = 1/4.001,muH = 0.15,muU = 0.4,xi = 0.53,omegaU = 0.29,omegaH=0.14,E0=None,Ia0=None,Is0=None,H0=None,U0=None,R0=None,D0=None, fittingByCumulativeCases=True, yWeight=1,dWeight = 1,hosWeight=1,uWeight=1,stand_error = True, particles = 300, itera = 1000, c1 = 0.1, c2 = 0.3, w = 0.9, k = 5, norm = 2,standarPSO=True):
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
        
        if E0==None:
            E0=[0,y[0]*10/self.N]
        if Ia0==None:
            Ia0=[0,y[0]*10/self.N]
        if Is0==None:
           Is0=[0,y[0]/self.N]
        if H0==None:
            if hos:
                H0=hos[0]
            else:
                H0 = 0
        if U0==None:
            if u:
                U0 = u[0]
            else:
                U0 = 0
        if R0==None:
            R0 = 0
        if D0 == None:
            D0 = d[0]
        
        self.coef = self.Coefficient(BetaChange,standarPSO)
        self.BetaChange = BetaChange
        if not self.coef.addBeta(beta):
                return
        if self.BetaChange!=0:
            if dayBetaChange==None:
                dayBetaChange = []
                for i in range(BetaChange):
                    dayBetaChange.append([])
                    dayBetaChange[i].append(x[5])
                    dayBetaChange[i].append(x[-4])
                if not self.coef.addDay(dayBetaChange):
                    return
            else:
                if not self.coef.addDay(dayBetaChange):
                    return
            if not self.coef.addDay(dayBetaChange):
                return

        if not self.coef.addVar('gammaH',gammaH):
            return
        if not self.coef.addVar('gammaU',gammaU):
            return
        if not self.coef.addVar('delta',delta):
            return
        if not self.coef.addVar('kappa',kappa):
            return
        if not self.coef.addVar('h',h):
            return
        if not self.coef.addVar('p',p):
            return
        if not self.coef.addVar('gammaA',gammaA):
            return
        if not self.coef.addVar('gammaS',gammaS):
            return
        if not self.coef.addVar('muH',muH):
            return
        if not self.coef.addVar('muU',muU):
            return
        if not self.coef.addVar('xi',xi):
            return
        if not self.coef.addVar('omegaU',omegaU):
            return
        if not self.coef.addVar('omegaH',omegaH):
            return
        if not self.coef.addVar('E0',E0):
            return
        if not self.coef.addVar('Ia0',Ia0):
            return
        if not self.coef.addVar('Is0',Is0):
            return
        if not self.coef.addVar('H0',H0):
            return
        if not self.coef.addVar('U0',U0):
            return
        if not self.coef.addVar('R0',R0):
            return
        if not self.coef.addVar('D0',D0):
            return
#standarPSO=True
        self.y = np.array(y)
        self.d = np.array(d)
        self.x = np.array(x)
        if hos:
            self.hos = np.array(hos)
        else:
            self.hos = hos
        if u:
            self.u = np.array(u)
        else:
           self.u = u 
        self.fittingByCumulativeCases = fittingByCumulativeCases
        self.yWeight = yWeight
        self.dWeight = dWeight
        self.hosWeight = hosWeight
        self.uWeight = uWeight
        self.stand_error = stand_error
        self.particles = particles
        self.itera = itera
        self.c1 = c1
        self.c2 = c2
        self.w = w
        self.k = k
        self.norm = norm
        self.NC = self._Models__changeCases(self.y)
        self.standarPSO = standarPSO
        options = {'c1': c1, 'c2': c2, 'w': w,'k':k,'p':norm}
        dim = self.coef.getDimention()

        if dim==0:
            self.beta = beta
            self.dayBetaChange = dayBetaChange
            self.gammaH = gammaH
            self.gammaU = gammaU
            self.delta = delta
            self.kappa = kappa
            self.h = h
            self.p=p
            self.gammaA = gammaA
            self.gammaS = gammaS
            self.muH = muH
            self.muU = muU
            self.xi = xi
            self.omegaU = omegaU
            self.omegaH = omegaH
            self.E0 = E0
            self.Ia0 = Ia0
            self.Is0 = Is0
            self.H0 = H0
            self.U0 = U0
            self.R0 = R0
            self.D0 = D0
        
        else:
            optimizer = ps.single.LocalBestPSO(n_particles=particles, dimensions=dim, options=options,bounds=self.coef.getBound())
            cost, pos = optimizer.optimize(self._objectiveFunction, itera,n_processes=self.nCores)
            self.pos = pos
            dic = self.coef.getDic(pos)
            self.beta = dic['beta']
            if self.BetaChange>0:
                self.dayBetaChange = dic['dayBetaChange']
            self.gammaH = dic['gammaH']
            self.gammaU = dic['gammaU']
            self.delta = dic['delta']
            self.kappa = dic['kappa']
            self.h = dic['h']
            self.p=dic['p']
            self.gammaA = dic['gammaA']
            self.gammaS = dic['gammaS']
            self.muH = dic['muH']
            self.muU = dic['muU']
            self.xi = dic['xi']
            self.omegaU = dic['omegaU']
            self.omegaH = dic['omegaH']
            self.E0 = dic['E0']
            self.Ia0 = dic['Ia0']
            self.Is0 = dic['Is0']
            self.H0 = dic['H0']
            self.U0 = dic['U0']
            self.R0 = dic['R0']
            self.D0 = dic['D0']
            self.cost_history = optimizer.cost_history
        self.isFit=True
        self.predict(0)

    def predict(self,numDays):
        ''' x = dias passados do dia inicial 1'''
        if numDays<0:
            print('\nnumDays must be a positive number!\n')
            return
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        x = range(1,len(self.y)+1+numDays)
        self.xpred = x
        self.predictNumDays = numDays
        
        if self.BetaChange>0:
            S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO_2(x,self.dayBetaChange,self.beta,self.gammaH,self.gammaU,self.delta,self.kappa,self.h,self.p,self.gammaA,self.gammaS,self.muH,self.muU,self.xi,self.omegaU,self.omegaH,self.E0,self.Ia0,self.Is0,self.H0,self.U0,self.R0,self.D0)
        else:
            S,E,IA,IS,H,U,R,D,Nw = self.__cal_EDO(x,self.beta,self.gammaH,self.gammaU,self.delta,self.kappa,self.h,self.p,self.gammaA,self.gammaS,self.muH,self.muU,self.xi,self.omegaU,self.omegaH,self.E0,self.Ia0,self.Is0,self.H0,self.U0,self.R0,self.D0)
        self.ypred = Nw
        self.dpred = D
        self.S = S
        self.E = E
        self.IA = IA
        self.IS = IS
        self.H = H
        self.U = U
        self.R = R  
        self.NCpred =self._Models__changeCases(self.ypred)
        return self.ypred

    def plot(self,local=None,InitialDate=None,CompartmentPlots=None,SaveFile=None):
        
        if InitialDate != None:
            initial_date=date(int(InitialDate[0:4]), int(InitialDate[5:7]), int(InitialDate[8:11]))
            dates = []
            dates.append(initial_date.strftime('%Y-%m-%d'))  
            for i in range(len(self.ypred)-1):
                d=initial_date + timedelta(days=i)
                dates.append(d.strftime('%Y-%m-%d')) 
        else:
            dates=self.xpred
       
        
        #Plotting
        
        if CompartmentPlots==None:           
            fig, ax = plt.subplots(figsize=(17,10))
            ax.grid(which='major', axis='both', color='black',linewidth=1.,alpha=0.3)
            ax.plot(dates,self.ypred,'b-', linewidth=2.5,label='Model')
            ax.scatter(dates[:len( self.y)], self.y,  s=18,color='black',label='Reported data',zorder=3)
            if self.BetaChange >0:
                for day in self.dayBetaChange:
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')
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
            if self.BetaChange >0:
                for day in self.dayBetaChange:
                    ax.axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')
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
                gs = gridspec.GridSpec(nrows = 1, ncols = 4)
                fig=plt.figure(figsize=(2*15.4,10))
                gs.update(wspace = 0.55)
                gs.update(hspace = 0.55)
                ax=[]
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

                if self.BetaChange >0:
                    for day in self.dayBetaChange:
                        ax[k].axvline(day, 0, 600,c='r',linestyle='--',label='Beta change')

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
        

    def getCoef(self):
        res = {}
        if self.isFit:
            if self.BetaChange==0:
                res['beta']=self.beta[0]
            else:
                for i in range(len(self.beta)):
                    res['beta'+str(i)] = self.beta[i]
                for i in range(len(self.dayBetaChange)):
                    res['dayBetaChange'+str(i)] = self.dayBetaChange[i]
            #R0=None,D0=None
            res['gammaH'] = self.gammaH
            res['gammaU'] = self.gammaU
            res['delta'] = self.delta
            res['kappa'] = self.kappa
            res['h'] = self.h
            res['p'] = self.p
            res['gammaA'] = self.gammaA
            res['gammaS'] = self.gammaS
            res['muH'] = self.muH
            res['muU'] = self.muU
            res['xi'] = self.xi
            res['omegaU'] = self.omegaU
            res['omegaH'] = self.omegaH
            res['E0'] = self.E0
            res['Ia0'] = self.Ia0
            res['Is0'] = self.Is0
            res['H0'] = self.H0
            res['U0'] = self.U0
            res['R0'] = self.R0
            res['D0'] = self.D0
        else:
            print('The model is not fitted!\n')
            res = None
        return res
    
    def getEstimation(self):
        return dict(zip(['S','E','IA','IS','H','U','R','D','Cumulative_cases_predict','new_cases_predict'],[self.S,self.E,self.IA,self.IS,self.H,self.U,self.R,self.D,self.ypred,self.NCpred]))
    
    def computeCI(self, times=500, level=0.95):
        if self.isFit==False:
            print('\nModels is not fitted\n')
            return None
        if self.coef.getDimention()==0:
            print('none of coefficients are fitted\n')
            return
        if self.isCI:
            self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
            self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
            self.lE = self._Models__getConfidenceInterval(self._bE, level,False)
            self.lIA = self._Models__getConfidenceInterval(self._bIA, level,False)
            self.lIS = self._Models__getConfidenceInterval(self._bIS, level,False)
            self.lH = self._Models__getConfidenceInterval(self._bH, level,False)
            self.lU = self._Models__getConfidenceInterval(self._bU, level,False)
            self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
            self.lD = self._Models__getConfidenceInterval(self._bD, level,False)
            self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
            self.lCoef = {}
            for k in self.bcoef.keys():
                self.lCoef[k]=self._Models__getConfidenceInterval(self._bCoef[k], level,True)
                        
        #Define empty lists to recive results
        self._bypred = []
        self._bS = []
        self._bE = []
        self._bIA = []
        self._bIS = []
        self._bH = []
        self._bU = []
        self._bR = []
        self._bD = []
        self._bNCpred = []
        
        listCoef = tuple(self.getCoef().keys())
        self._bCoef={}
        for k in listCoef:
            self._bCoef[k]=[]
        casesSeries = self._Models__genBoot(self.y, times)
        copia = copy.deepcopy(self)
        for i in range(0,len(casesSeries)):
            copia.__fit(x=self.x, y=casesSeries[i])
            copiaCoef = copia.getCoef()
            for k in listCoef:
                self._bCoef[k].append(copiaCoef[k])
            self._bypred.append(copia.ypred)
            self._bS.append(copia.S)
            self._bE.append(copia.E)
            self._bIA.append(copia.IA)
            self._bIS.append(copia.IS)
            self._bH.append(copia.H)
            self._bU.append(copia.U)
            self._bR.append(copia.R)
            self._bD.append(copia.D)
            self._bNCpred.append(copia.NCpred)
                       
        
        self.lypred = self._Models__getConfidenceInterval(self._bypred, level,False)
        self.lS = self._Models__getConfidenceInterval(self._bS, level,False)
        self.lE = self._Models__getConfidenceInterval(self._bE, level,False)
        self.lIA = self._Models__getConfidenceInterval(self._bIA, level,False)
        self.lIS = self._Models__getConfidenceInterval(self._bIS, level,False)
        self.lH = self._Models__getConfidenceInterval(self._bH, level,False)
        self.lU = self._Models__getConfidenceInterval(self._bU, level,False)
        self.lR = self._Models__getConfidenceInterval(self._bR, level,False)
        self.lD = self._Models__getConfidenceInterval(self._bD, level,False)
        self.lNCpred = self._Models__getConfidenceInterval(self._bNCpred, level,False)
        self.lCoef = {}
        
        for k in listCoef:
            self.lCoef[k] = self._Models__getConfidenceInterval(self._bCoef[k], level,True)
        self.isCI=True        