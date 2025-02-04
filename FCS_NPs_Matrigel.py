# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 09:42:48 2024

@author: Alberto
"""



import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
from scipy.stats import ttest_rel,ttest_ind, shapiro, f_oneway, tukey_hsd,chisquare
from iapws import IAPWS97 
from scipy.optimize import curve_fit

from sklearn.metrics import r2_score 
from scipy.stats import linregress


os.chdir(r"XXXXXXXX")
    
    



sns.set_theme()
sns.set_style("ticks")
sns.set_context("paper",font_scale=1.5)
sns.set_palette("bright")



def diffusion(tau,g0,td):
    

    y = g0/((1+(tau/td))*(1+(w**2/z0**2)*(tau/td))**(1/2))
    
    return y

def anomalous_D(tau,g0,td,a,):
    
    y = g0/((1+(tau/td)**a)*(1+(w**2/z0**2)*(tau/td)**a)**(1/2))
    
    return y
    

kb = 1.38e-23


T = 273.15+36

#Calibrations for every repetition

ws = {"100424":0.24816,"160424":0.24862,"140324":0.24862,"170424":0.23252}

z0s = {"100424":8.68965,"160424":4.60258,"140324":4.60258,"170424":2.18211}

D_values = pd.DataFrame()
D_pure_values = pd.DataFrame()

together = pd.DataFrame()

k = 0
n = 0
j = 0

#Data analysis
for folder in os.listdir():
    
    if "Matrigel" in folder:
        
        sample = str(folder.split(" ")[1]) +" mg/ml"
        
        rep = 1
        
        for day in os.listdir(f"{folder}"):
            if day == "100424" and sample != "3.0 mg/ml":
                pass 
            else:
                chip = 1
                w = ws[day]
                z0 = z0s[day]
                for file in os.listdir(f"{folder}/{day}"):
                    
                    data = pd.read_csv(f"{folder}/{day}/{file}",skiprows=1,sep="\t",skipinitialspace=True)
                    
                    for col in data.columns:
                        if "Correlation Channel" in col:
                            
                            
                                G_exp = data.loc[:,col].dropna()
                                
                                t = data.loc[:len(G_exp)-1,"Time [ms]"]
                                
                                pars,cov = curve_fit(anomalous_D,t,G_exp)
                                
                                (g0,tau_d,a) = pars
                                
                                d_w = 2/a
                                
                                D = w**2/(4*tau_d*1e-3)
                                
                                G_fit = anomalous_D(t,*pars)
                                
                                chi_2 = r2_score(G_exp,G_fit)
                                
                                
                                pars_pure,cov_pure = curve_fit(diffusion,t,G_exp)
                                
                                (g0_pure,tau_d_pure) = pars_pure
                                D_pure = w**2/(4*tau_d_pure*1e-3)
                                
                                if chi_2>0.99:
                                    
                                    
                                    
                                    
                                    D_values.loc[k,"Matrigel Concentration"] = sample
                                    D_values.loc[k,"Repetition"] = day
                                    D_values.loc[k,"Chip"] = chip
                                    D_values.loc[k,"D ($\mu$m$^2$/s)"] = D
                                    D_values.loc[k,"G$_0$"] = g0
                                    D_values.loc[k,"Fractal dimension of diffusion (d$_w$)"] = 2/a
                                    D_values.loc[k,r"Diffusion Exponent ($\alpha$)"] = a
                                    
                                    together.loc[j,"Matrigel Concentration"] = sample
                                    together.loc[j,"D ($\mu$m$^2$/s)"] = D
                                    together.loc[j,"G$_0$"] = g0
                                    together.loc[j,r"Diffusion Exponent ($\alpha$)"] = a
                                    together.loc[j,r"Model"] = "Anomalous"
                                    together.loc[j,"Repetition"] = day
                                    together.loc[j,"Chip"] = chip
                                    
                                    j+=1
                                    k+=1
                                
                                G_fit_pure = diffusion(t, *pars_pure)
                                
                                
                                
                                D_pure_values.loc[k,"Matrigel Concentration"] = sample
                                D_pure_values.loc[k,"D$_{pure}$ ($\mu$m$^2$/s)"] = D_pure
                                D_pure_values.loc[k,"G$_0$ pure"] = g0_pure
                                D_pure_values.loc[k,"Repetition"] = day
                                D_pure_values.loc[k,"Chip"] = chip
                                
                                together.loc[j,"Matrigel Concentration"] = sample
                                together.loc[j,"D ($\mu$m$^2$/s)"] = D_pure
                                together.loc[j,"G$_0$"] = g0_pure
                                together.loc[j,r"Diffusion Exponent ($\alpha$)"] = 1
                                together.loc[j,r"Model"] = "Free"
                                together.loc[j,"Repetition"] = day
                                together.loc[j,"Chip"] = chip
                                
                                
                                
                                j+=1
                                
                                
                                if n == 164:
                                    
                                    
                                    plt.plot(t,G_exp,".",label="ACF")
                                    plt.plot(t,G_fit,"-",color = "#e8000b",label = "Anomalous Diffusion")
                                    plt.plot(t,G_fit_pure,"-",color="#1ac938",label="Free Diffusion")
                                    #plt.title(f"{day} {sample} {rep} {chip}")
                                    plt.xscale("log")
                                    plt.xlabel("Time (ms)")
                                    plt.ylabel(r"G ($\tau$)")
                                    plt.legend()
                                    os.makedirs(f"Figures/G_curves/{folder}/{day}",exist_ok=True)
                                    
                                    plt.savefig(f"Figures/G_curves/{folder}/{day}/{file}_{n}.tif",dpi=100,bbox_inches="tight")
                                    plt.show()
                                    plt.close()
                                n+=1
                                
                    
                    chip+=1
                rep+=1                

#Plotting and statistics
means = D_values.groupby(by = ["Matrigel Concentration","Repetition","Chip"],as_index=False).mean()
print(means)

means2 = means.groupby(by = ["Matrigel Concentration","Repetition"],as_index=False).mean()


sns.barplot(data=means2,x = "Matrigel Concentration",y = "D ($\mu$m$^2$/s)",palette = "tab10",errorbar="se",capsize=0.1,saturation=1)
plt.ylabel(r"Diffusion Coefficient ($\mu$m$^2$/s$^{\alpha}$)")
plt.savefig("Figures/Diffusion Coefficients.tif",bbox_inches="tight",dpi=300)

plt.show()
plt.close()

sns.histplot(data = D_values,x = "D ($\mu$m$^2$/s)",hue = "Matrigel Concentration",bins=20,kde=True,palette = "tab10")
plt.savefig("Figures/Diffusion Coefficients histograms.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()

sns.barplot(data=means2,x = "Matrigel Concentration",y = "Fractal dimension of diffusion (d$_w$)",palette = "tab10",errorbar="se",capsize=0.1,saturation=1)

plt.savefig("Figures/Fractal Dimension of Diffusion.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()

sns.barplot(data=means2,x = "Matrigel Concentration",y = r"Diffusion Exponent ($\alpha$)",palette = "hls",errorbar="se",capsize=0.1,saturation=1)

plt.savefig("Figures/Diffusion exponent.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()


sns.barplot(data=means2,x = "Matrigel Concentration",y = "G$_0$",palette = "Set2",errorbar="se",capsize=0.1,saturation=1)

plt.savefig("Figures/g0 anomalous.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()

Ds = []
Ds_err = []
g0s = []
g0s_err = []

for media in set(means2.loc[:,"Matrigel Concentration"]):
    
    Ds.append(np.mean(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"D ($\mu$m$^2$/s)"]))
    Ds_err.append(np.std(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"D ($\mu$m$^2$/s)"])/np.sqrt(len(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"D ($\mu$m$^2$/s)"])))
    g0s.append(np.mean(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"G$_0$"]))
    g0s_err.append(np.std(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"G$_0$"])/np.sqrt(len(means2.loc[means2.loc[:,"Matrigel Concentration"]==media,"D ($\mu$m$^2$/s)"])))
  

plt.errorbar(x = Ds,y = g0s,xerr = Ds_err,yerr = g0s_err,fmt="k.",barsabove=True,capsize=2.5,label="Experimental")
plt.xlabel("D ($\mu$m$^2$/s)")
plt.ylabel("G$_0$")

regression_result = linregress(np.log10(Ds),np.log10(g0s))

m = regression_result.slope
b = regression_result.intercept
rvalue = regression_result.rvalue
pvalue = regression_result.pvalue
m_err = regression_result.stderr
b_err = regression_result.intercept_stderr

plt.plot(Ds,10**(m*np.log10(Ds)+b),"r",label="Fit")

plt.legend()
plt.savefig("Figures/g0 vs D.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()
    
means_pure = D_pure_values.groupby(by = ["Matrigel Concentration","Repetition","Chip"],as_index=False).mean()
print(means)

means2_pure = means_pure.groupby(by = ["Matrigel Concentration","Repetition"],as_index=False).mean()


sns.barplot(data=means2_pure,x = "Matrigel Concentration",y = "D$_{pure}$ ($\mu$m$^2$/s)",palette = "tab10",errorbar="se",capsize=0.1,saturation=1)
plt.ylabel("Diffusion Coefficient (Pure) ($\mu$m$^2$/s)")
plt.savefig("Figures/Diffusion Coefficients Pure.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()

sns.barplot(data=means2_pure,x = "Matrigel Concentration",y = "G$_0$ pure",palette = "Set2",errorbar="se",capsize=0.1,saturation=1)

plt.savefig("Figures/g0 Pure.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()



sns.histplot(data = D_pure_values,x = "D$_{pure}$ ($\mu$m$^2$/s)",hue = "Matrigel Concentration",bins=20,kde=True,palette = "tab10")
plt.savefig("Figures/Pure Diffusion Coefficients histogram.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()



populations_anomalous = []

names = []
for media in set(means2.loc[:,"Matrigel Concentration"]):
    
    
    condition = list(means2.loc[means2.loc[:,"Matrigel Concentration"] == media,"D ($\mu$m$^2$/s)"])
    
    names.append(media)
    
    populations_anomalous.append(condition)

print(names)
print(populations_anomalous)

stats_anomalous = tukey_hsd(*populations_anomalous)

print(stats_anomalous)



populations_pure = []

names = []
for media in set(means2_pure.loc[:,"Matrigel Concentration"]):
    
    
    condition = list(means2_pure.loc[means2_pure.loc[:,"Matrigel Concentration"] == media,"D$_{pure}$ ($\mu$m$^2$/s)"])
    
    names.append(media)
    
    populations_pure.append(condition)

print(names)
print(populations_pure)

stats_pure = tukey_hsd(*populations_pure)

print(stats_pure)



populations_fractal = []

names = []
for media in set(means2.loc[:,"Matrigel Concentration"]):
    
    
    condition = list(means2.loc[means2.loc[:,"Matrigel Concentration"] == media,r"Diffusion Exponent ($\alpha$)"])
    
    names.append(media)
    
    populations_fractal.append(condition)

print(names)
print(populations_fractal)

stats_fractal = tukey_hsd(*populations_fractal)

print(stats_fractal)



G0s = []
G0s_err = []
Ds = []
Ds_err = []
alphas = []
alphas_err = []

Cs = list(set(means2.loc[:,"Matrigel Concentration"]))

Cs.sort()

Cs_n = [float(C.split(" ")[0]) for C in Cs]

for C in Cs:
    G0s.append(np.mean(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,"G$_0$"])))
    G0s_err.append(np.std(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,"G$_0$"])))
    Ds.append(np.mean(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,"D ($\mu$m$^2$/s)"])))
    Ds_err.append(np.std(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,"D ($\mu$m$^2$/s)"]))/np.sqrt(len(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,"D ($\mu$m$^2$/s)"]))))
    alphas.append(np.mean(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,r"Diffusion Exponent ($\alpha$)"])))
    alphas_err.append(np.std(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,r"Diffusion Exponent ($\alpha$)"]))/np.sqrt(len(list(means2.loc[means2.loc[:,"Matrigel Concentration"]==C,r"Diffusion Exponent ($\alpha$)"]))))
    


plt.errorbar(Cs_n,G0s,yerr = G0s_err,fmt="k.-",barsabove=True,capsize=2.5)
plt.ylabel("G$_0$")
plt.xlabel("Matrigel Concentration (mg/ml)")
plt.show()
plt.close()


plt.errorbar(Cs_n,Ds,yerr = Ds_err,fmt="b.-",barsabove=True,capsize=2.5)
plt.ylabel("Diffusion Coefficient ($\mu$m$^2$/s)")
plt.xlabel("Matrigel Concentration (mg/ml)")

ax = plt.gca()

ax2 = ax.twinx()


ax2.errorbar(Cs_n,alphas,yerr = alphas_err,fmt="r.-",barsabove=True,capsize=2.5)
ax2.set_ylabel(r"Diffusion Exponent ($\alpha$)")

ax.yaxis.label.set_color('blue')
ax.tick_params(axis='y', colors='blue')

ax2.yaxis.label.set_color('red')
ax2.tick_params(axis='y', colors='red')


plt.show()
plt.close()






means_together = together.groupby(by = ["Matrigel Concentration","Repetition","Model"],as_index=False).mean()

sns.barplot(data = means_together,x = "Matrigel Concentration", y = "D ($\mu$m$^2$/s)", hue = "Model",palette = "tab10",errorbar = "se",saturation=1,capsize=0.1)
plt.show()
plt.close()





ROI_means = means2.loc[:,["Matrigel Concentration","D ($\mu$m$^2$/s)","G$_0$",r"Diffusion Exponent ($\alpha$)"]]








