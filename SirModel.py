from scipy.integrate import odeint
from scipy.signal import StateSpace, lsim
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import mpld3


#import casadi  as csdi
#x = csdi.MX.sym("x")
#print(csdi.jacobian(csdi.sin(x),x))

class FHMData:
    def __init__(self):
        # https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/bekraftade-fall-i-sverige/
        filename = "Folkhalsomyndigheten_Covid19.xlsx"
        self.data = pd.read_excel(filename,sheet_name="Antal avlidna per dag", index_col=[0])[:-1]

class SirModel:
#
#    S - I - R model
#
#    Susceptibles   S(t)
#    Infected       I(t)
#    Recovered      R(t)  
#    Death          D(t)
#   
#    Population N(t) = S(t) +  I(t) + R(t) - D(t)
#    
#    Initial condition
#    N(0)  =  No  = 10 000 001
#    S(0)  =  So  = 10 000 000
#    I(0)  =  Io  = 1
#    R(0)  =  Ro  = 0
#
#    dS/dt = -k01 * I * S                       (0)
#    dI/dt =  k01 * I * S - k12 * I - k13 * I   (1)
#    dR/dt =  k12 * I                           (2)
#    dD/dt =  k13 * I                           (3)
#
#    Ro    = k01 * So / k12
#
#    r      = 2^(1/3)
#    Re     = 2.5
#    gamma  = (r - 1) / (Re - 1) 
#    lambda = r - 1 + gamma
#    k01    = lambda / N
#    k12    = gamma

#    https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

#    https://www.medrxiv.org/content/10.1101/2020.04.15.20066050v1.full.pdf

    def __init__(self, Ro = 2.5, k12=0.325, k13=0.00011, So = 1E7, dateStart = '2020-01-22', plotDateRange = ['2020-03-01','2020-06-01']):
        self.startDate = datetime.datetime.strptime(dateStart, '%Y-%m-%d')  
        self.plotStartDate = datetime.datetime.strptime(plotDateRange[0], '%Y-%m-%d') #datetime.date(2020,3,1)
        self.plotEndDate   = datetime.datetime.strptime(plotDateRange[1], '%Y-%m-%d') #datetime.date(2020,6,1)
        self.So        = So
        self.Io        = 1
        self.k12       = self.interp(k12)
        self.k01       = self.calcK01(Ro)
        #self.k01       = self.interp(k01 / So)        # k01=0.52,
        
        self.k13       = self.interp(k13)

        self.Ro        = self.k01(0) * So / k12
        self.dti       = pd.date_range(dateStart, periods=(self.plotEndDate - self.startDate).days, freq='D')
        self.t         =  (self.dti-self.dti.min()).astype('timedelta64[D]').astype(int)
        self.simResult     = self.solve()

        #plt.plot(self.k01(self.dti) / k12(self.dti) * self.simResult['Susceptibles'])
        #plt.show()
        self.plot(self.t, self.simResult)
        for key in Ro:
            print(f"Data {key} Ro: {Ro[key]:5.2f}")

    @staticmethod
    def interp(k, t = [0 ,1]):
        return interpolate.interp1d(t,[k,k], bounds_error=False, fill_value=(k,k))

    def calcK01(self,Ro):
        if isinstance(Ro,dict):
            t , y = [], []
            for key in Ro:
                date = datetime.datetime.strptime(key, '%Y-%m-%d')
                delta = date.date() - self.startDate.date()
                t.append(delta.days)
                y.append((Ro[key]) / self.So * self.k12(delta.days)) 
            return interpolate.interp1d(t,y, bounds_error=False, fill_value=(y[0],y[-1]),kind='previous')
        else:
            k = Ro / self.So * self.k12(0)
            return interpolate.interp1d([0 , 1],[k,k], bounds_error=False, fill_value=(k,k))

            

    @staticmethod
    def f(x, t, k01, k12, k13):
        dx_dt = [0, 0, 0, 0]
        dx_dt[0] = -k01(t) * x[1] * x[0] 
        dx_dt[1] =  k01(t) * x[1] * x[0] - k12(t) * x[1] - k13(t) * x[1] 
        dx_dt[2] =  k12(t) * x[1]
        dx_dt[3] =  k13(t) * x[1] 
        return dx_dt

    def solve(self):
        s  = odeint(self.f, y0=[self.So, self.Io, 0, 0], t=self.t, args=(self.k01, self.k12, self.k13))
        simResult = {}
        simResult['Susceptibles'] = pd.Series(s[:,0], self.dti)
        simResult['Infected']     = pd.Series(s[:,1], self.dti)
        simResult['Recovered']    = pd.Series(s[:,2], self.dti)
        simResult['Death']        = pd.Series(s[:,3], self.dti)
        return simResult

    def plot(self,t, simResult):
        plt.figure(0,figsize=(15,15))
        x, y = (4, 3)
        plt.subplot(x,y,1)
        plt.plot( simResult['Susceptibles'] )
        plt.title(f'Friska (Ro: {self.Ro:5.2f})')
        plt.ylabel('Antal')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,2)
        plt.plot( simResult['Susceptibles'] )
        plt.yscale('log')
        plt.title('Friska')
        plt.ylabel('Antal')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)
        
        plt.subplot(x,y,3)
        plt.plot( simResult['Susceptibles'].diff() )
        plt.title('Friska')
        plt.ylabel('Antal per dag')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,4)
        plt.plot(simResult['Infected'])
        plt.title('Sjuka')
        plt.ylabel('Antal')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)
        
        plt.subplot(x,y,5)
        plt.plot(simResult['Infected'])
        plt.yscale('log')
        plt.title('Sjuka')
        plt.ylabel('Antal')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.ylim((1))
        plt.grid(True)

        plt.subplot(x,y,6)
        plt.plot(simResult['Infected'].diff())
        plt.title('Sjuka')
        plt.ylabel('Antal per dag')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,7)
        plt.plot(simResult['Recovered'])
        plt.title('Friska immuna')
        plt.ylabel('Antal')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,8)
        plt.plot(simResult['Recovered'])
        plt.yscale('log')
        plt.title('Friska immuna')
        plt.ylabel('Antal')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.ylim((1))
        plt.grid(True)
        
        plt.subplot(x,y,9)
        plt.plot(simResult['Recovered'].diff())
        plt.title('Friska immuna')
        plt.ylabel('Antal per dag')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        fhm = FHMData()
        plt.subplot(x,y,10)
        plt.plot(simResult['Death'])
        plt.plot(fhm.data.Antal_avlidna.cumsum())
        plt.title('Döda')
        plt.xlabel('Dagar')
        plt.ylabel('Antal')
        plt.legend(['Sim','FHM'])
        plt.xlim([self.plotStartDate ,self.plotEndDate])
        plt.gcf().autofmt_xdate()
        plt.grid(True)

        plt.subplot(x,y,11)
        plt.plot(simResult['Death'])
        plt.plot(fhm.data.Antal_avlidna.cumsum())
        plt.yscale('log')
        plt.title('Döda')
        plt.xlabel('Dagar')
        plt.ylabel('Antal')
        plt.ylim(1)
        plt.gcf().autofmt_xdate()
        plt.grid(True)

        plt.subplot(x,y,12)
        plt.plot(simResult['Death'].diff())
        plt.plot(fhm.data.Antal_avlidna)
        plt.title('Döda')
        plt.xlabel('Dagar')
        plt.ylabel('Antal per dag')
        plt.xlim([self.plotStartDate ,self.plotEndDate])
        plt.grid(True)
        plt.gcf().autofmt_xdate()
        plt.tight_layout(pad=0.2, w_pad=0.1, h_pad=0.1)
        plt.subplots_adjust(wspace=0.25, hspace=0.15)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        #plt.savefig('covid-19_sim.png', bbox_inches='tight')
        #mpld3.show()
        plt.show()


Ro = {'2020-01-01': 2.4 ,'2020-03-15': 1.6 , '2020-04-03': 1.13, '2020-08-01':1.7}
sirm = SirModel(Ro = Ro, k12=0.325*0.95, k13=0.00023*2.2, So = 10E6, dateStart = '2020-02-24', plotDateRange = ['2020-03-01','2021-02-01'])
