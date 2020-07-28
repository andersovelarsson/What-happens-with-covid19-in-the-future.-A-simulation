from scipy.integrate import odeint
from scipy.signal import StateSpace, lsim
from scipy import interpolate
from scipy.optimize import minimize

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import mpld3
import urllib.request
import os, time



#import casadi  as csdi
#x = csdi.MX.sym("x")
#print(csdi.jacobian(csdi.sin(x),x))

class FHMData:
    def __init__(self):
        filename = "Folkhalsomyndigheten_Covid19.xlsx"
        modified = os.path.getmtime(filename)
        # https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/bekraftade-fall-i-sverige/
        if modified < time.time()-60*60*24:
            urllib.request.urlretrieve ("https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data", filename)

        self.data = pd.read_excel(filename,sheet_name="Antal avlidna per dag", index_col=[0])[:-1]

class SirModel:
#
#    S - I - R - D - model
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

    def __init__(self, Ro = {'2020-02-01':2.5}, k12=0.325, k13=0.00011, So = 1E7, dateStart = '2020-01-22', plotDateRange = ['2020-03-01','2020-06-01']):
        self.FHMData       = FHMData()
        self.startDate     = datetime.datetime.fromisoformat(dateStart)  
        self.plotStartDate = datetime.datetime.fromisoformat(plotDateRange[0])
        self.plotEndDate   = datetime.datetime.fromisoformat(plotDateRange[1])
        self.So            = So
        self.Io            = 1
        self.No            = self.So + self.Io                     # Initial population
        self.k12           = self.interp(k12)       
        self.k13           = self.interp(k13)
        self.dti           = pd.date_range(dateStart, periods=(self.plotEndDate - self.startDate).days, freq='D')
        self.t             = (self.dti-self.dti.min()).astype('timedelta64[D]').astype(int)
        self__Ro           = None
        self.Ro            = Ro
        
    @property
    def Ro(self):
        return self.__Ro

    @Ro.setter
    def Ro(self, Ro):
        def calcRo(Ro,startDate):
            if isinstance(Ro, interpolate.interpolate.interp1d):
                return Ro
            if isinstance(Ro,dict):
                t , y = [], []
                for key in Ro:
                    date  = datetime.datetime.fromisoformat(key)
                    delta = date.date() - startDate.date()
                    t.append(delta.days)
                    y.append(Ro[key]) 
                if len(Ro) == 1:
                    t.append(100)
                    y.append(list(Ro.values())[0])

                return interpolate.interp1d(t,y, bounds_error=False, fill_value=(y[0],y[-1]),kind='previous')
            else:
                k = Ro
                return interpolate.interp1d([0 , 1],[k,k], bounds_error=False, fill_value=(k,k))  
        
        def RoTabel(Ro):
            t , y = [], []
            s = pd.Series(self.__Ro.y,self.__Ro.x)
            for index, value in s.items():
                date = (self.startDate + datetime.timedelta(days=index)).strftime('%Y-%m-%d %H:%M:%S')
                #print(f"Date {date} Ro: {value:5.2f}")
                t.append(date)
                y.append(value)
            t.append(self.plotEndDate)
            y.append(y[-1])
            return pd.Series(y,t).to_string()

        self.__Ro = calcRo(Ro,self.startDate)
        self.RoText = RoTabel(Ro)
        
        self.simResult = self.solve(Ro = self.Ro, k12 = self.k12, k13 = self.k13, So = self.So)

    @staticmethod
    def interp(k, t = [0 ,1]):
        return interpolate.interp1d(t,[k,k], bounds_error=False, fill_value=(k,k))

    @staticmethod
    def sirdModel(x, t, Ro, k12, k13, So):
        dx_dt = [0, 0, 0, 0]
        R01 = Ro(t) * k12(t) / So * x[1] * x[0]
        R12 = k12(t) * x[1]
        R13 = k13(t) * x[1] 
        dx_dt[0] = -R01 
        dx_dt[1] =  R01 - R12 - R13 
        dx_dt[2] =  R12
        dx_dt[3] =  R13 
        return dx_dt

    def solve(self ,Ro ,k12 ,k13 ,So):
        s  = odeint(self.sirdModel, y0=[self.So, self.Io, 0, 0], t=self.t, args=(Ro, k12, k13, So))
        simResult = {}
        simResult['Susceptibles'] = pd.Series(s[:,0], self.dti)
        simResult['Infected']     = pd.Series(s[:,1], self.dti)
        simResult['Recovered']    = pd.Series(s[:,2], self.dti)
        simResult['Death']        = pd.Series(s[:,3], self.dti)
        simResult['Population']   = self.No - simResult['Death']
        return simResult

    def autoModelCalibration(self):
        def R(Ro ,x):
            Ro.y[0] = x[0]
            Ro.y[1] = x[1]
            Ro.y[2] = x[2]
            Ro.y[3] = x[3]
            Ro.y[4] = x[4]
            return Ro
        
        def residual(x):
            k12           = self.k12
            k13           = self.k13 
            measureData   = self.FHMData.data.Antal_avlidna.cumsum()
            Ro            = R(self.Ro, x)
            rmLastDays    = 10
            measureValues = measureData[:-rmLastDays]
            simResult     = self.solve(Ro = Ro ,k12 = k12, k13 = k13, So = self.So)
            simValues     = simResult['Death'].reindex(index=measureData.index[:-rmLastDays])
            return np.linalg.norm(measureValues - simValues)

        x = [2.4,1.6,1.11,1.2,1.35]
        res = minimize(residual, x, method='nelder-mead', options={'xatol': 1e-8, 'disp': True})
        print(res)
        return R(self.Ro ,res.x)

    def plot(self):
        t = self.t
        simResult = self.simResult
        plt.figure(figsize=(15,15))
        x, y = (4, 3)
        plt.subplot(x,y,1)
        plt.plot( simResult['Susceptibles'] / simResult['Population'] * 100.0 )
        plt.title(f'Susceptibles')
        plt.ylabel('Part of population [%]')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,2)
        plt.plot( simResult['Susceptibles'])
        plt.yscale('log')
        plt.title('Susceptibles')
        plt.ylabel('Number of people')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)
        
        ax = plt.subplot(x,y,3)
        #plt.step(self.Ro.index,self.Ro.values,where='post')
        plt.step(self.dti,self.Ro(self.t),where='post')
        plt.title(f'Reproduction number')
        plt.ylabel('Ro [-]')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gcf().autofmt_xdate()
        plt.grid(True)
        ax.text(0.95, 0.93, self.RoText,
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=10,bbox=dict(boxstyle="square",ec=(1., 1., 1.),fc=(1., 1., 1.), alpha=0.6) )
        
        plt.subplot(x,y,4)
        plt.plot(simResult['Infected'] / simResult['Population'] * 100.0 )
        plt.title('Infected')
        plt.ylabel('Part of population [%]')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)
        
        plt.subplot(x,y,5)
        plt.plot(simResult['Infected'])
        plt.yscale('log')
        plt.title('Infected')
        plt.ylabel('Number of people')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.ylim((1))
        plt.grid(True)

        #k01 * So / k12
        plt.subplot(x,y,6)
        plt.plot(self.dti,self.Ro(self.t)*simResult['Susceptibles'].values/self.So)
        plt.title(f'Effective reproduction number')
        plt.ylabel('R [-]')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gcf().autofmt_xdate()
        plt.grid(True)

        plt.subplot(x,y,7)
        plt.plot(simResult['Recovered'] / simResult['Population'] * 100.0)
        plt.title('Recovered')
        plt.ylabel('Part of population [%]')
        plt.xlim(plt.xlim([self.plotStartDate ,self.plotEndDate]))
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.grid(True)

        plt.subplot(x,y,8)
        plt.plot(simResult['Recovered'])
        plt.yscale('log')
        plt.title('Recovered (Sw: Friska immuna)')
        plt.ylabel('Number of people')
        plt.gca().axes.get_xaxis().set_major_formatter(plt.NullFormatter())
        plt.ylim((1))
        plt.grid(True)

        plt.subplot(x,y,10)
        plt.plot(simResult['Death'] / self.So * 1E6)
        plt.plot(self.FHMData.data.Antal_avlidna.cumsum() / self.So * 1E6)
        plt.title('Cumulative deaths per million people')
        plt.ylabel('Deaths per million')
        plt.legend(['Simulation','FHM Sweden'])
        plt.xlim([self.plotStartDate ,self.plotEndDate])
        plt.gcf().autofmt_xdate()
        plt.grid(True)

        plt.subplot(x,y,11)
        plt.plot(simResult['Death'])
        plt.plot(self.FHMData.data.Antal_avlidna.cumsum())
        plt.yscale('log')
        plt.title('Deaths')
        plt.ylabel('Number of deaths')
        plt.ylim(1)
        plt.gcf().autofmt_xdate()
        plt.grid(True)

        plt.subplot(x,y,12)
        plt.plot(simResult['Death'].diff())
        plt.plot(self.FHMData.data.Antal_avlidna)
        plt.title('Deaths')
        plt.ylabel('Daily deaths')
        plt.xlim([self.plotStartDate ,self.plotEndDate])
        plt.grid(True)
        plt.gcf().autofmt_xdate()
        plt.tight_layout(pad=0.2, w_pad=0.1, h_pad=0.1)
        plt.subplots_adjust(wspace=0.25, hspace=0.15)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        #plt.savefig('covid-19_sim.png', bbox_inches='tight')
        #mpld3.show()

        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()

        plt.show()
    def plot2(self):
        pass

if __name__ == "__main__":
    sirdm = SirModel(k12=0.3077, k13=0.000506, So = 10E6, dateStart = '2020-02-24', plotDateRange = ['2020-03-01','2022-03-01'])
    sirdm.Ro = {'2020-01-01': 2.37713217,
                '2020-03-16': 1.62439360, 
                '2020-04-02': 1.10813591, 
                '2020-04-24': 1.20101268, 
                '2020-05-23': 1.33650182,
                '2020-08-30': 1.33650182} 
    sirdm.plot()
    sirdm.Ro  = sirdm.autoModelCalibration()
    sirdm.plot()

