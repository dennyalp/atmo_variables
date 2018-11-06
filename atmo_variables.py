# -*- coding: utf-8 -*-
#"""
#Created on Thu Aug 16 15:07:07 2018
#
#@author: denny
#"""
# =============================================================================
# Atmospheric Science varaibles module. Mainly for ABL data analysis. 
#
# =============================================================================
import numpy as np
#==============================================================================
def pot_temp(pres,tair):
    
    """
    Estimates Potential temperature (K) 
    from Pressure (hPa) and Temperature (°C) 
    """
    
    tairk =  tair + 273.15 # converting to Kelvin
    theta =  tairk*((1000/pres)**0.286) 

    return theta

#==============================================================================
def sat_vap_pres(tair):
    
    """
    Estimates Saturated vapor pressure (hPa) 
    from Temperature (°C) 
    """
    tairk = tair+273.15
    E0    = 0.611      # kPa 
    B1    = 17.2694    # 1/K
    T1    = 273.16     # K
    T2    = 35.86      # K
    
    VAR1  = tairk - T1
    VAR2  = tairk - T2
    VAR3  = B1*VAR1
    satvp = E0 * np.exp(VAR3/VAR2)*10# TETEN'S FORMULA FOR SAT. VAPOUR PRES. (hPa)
    
    return satvp

#=============================================================================
def sp_hum(pres,tair,rh):
    """
    Estimates Specific humidity (g/kg) 
    from Pressure (hPa), Temperature (°C) and Relative humidity (%)
    """
    tairk = tair+273.15;
    rhh   = rh/100;
    pres  = pres/10;
    
    E0    = 0.611;      # kPa 
    B1    = 17.2694;    # 1/K
    T1    = 273.16;     # K
    T2    = 35.86;      # K
    
    VAR1  = tairk - T1;
    VAR2  = tairk - T2;
    VAR3  = B1*VAR1;
    satvp = E0 * np.exp(VAR3/VAR2);         # TETEN'S FORMULA FOR SAT. VAPOUR PRES. (kPa)
    vp    = satvp*rhh;                      # VAPOR PRES (kPa)   
    
    EPS   = 0.622;              # RD(-----)/RV(461 J/Kkg);  
    #VAR4  = EPS*satvp;
    #VAR5  = pres - satvp;
    #SATMR = VAR4/VAR5;          # SAT. MIXING RATIO IN g/g
    
    VAR6   = EPS*vp;
    SH_bfr = VAR6/pres;         # Specific humidity IN g/g
    SHUM   = SH_bfr*1000;         # Sp.hum in g/kg
    return SHUM

#==============================================================================

def mix_ratio(pres,tair,rh):
    """
    Estimates Mixing ratio (g/kg) 
    from Pressure (hPa), Temperature (°C) and Relative humidity (%)
    """
    tairk = tair+273.15;
    rhh   = rh/100;
    pres  = pres/10;
    
    E0    = 0.611;      # kPa 
    B1    = 17.2694;    # 1/K
    T1    = 273.16;     # K
    T2    = 35.86;      # K
    
    VAR1  = tairk - T1;
    VAR2  = tairk - T2;
    VAR3  = B1*VAR1;
    satvp = E0 * np.exp(VAR3/VAR2);      # TETEN'S FORMULA FOR SAT. VAPOUR PRES. (kPa)
    vp    = satvp*rhh;                   # VAPOR PRES (kPa)   
    
    EPS    = 0.622;              # RD(-----)/RV(461 J/Kkg);  
    VAR6   = EPS*vp;
    VAR7   = pres - vp
    MR_bfr = VAR6/VAR7      # MIXING RATIO IN g/g
    MR     = MR_bfr*1000;   # mix ratio in g/kg
    return MR
#==============================================================================
    
def vir_pot_temp(pres,tair,rh):
    
    """
    Estimates virtual potential temperature (K) 
    from Pressure (hPa), Temperature (°C) and Relative humidity (%)
    """
    
    MR=mix_ratio(pres,tair,rh) # mixing ratio in g/kg
    mr_bfr=MR/1000 #mixing ratio in g/g
    THETA=pot_temp(pres,tair);
    THETAV   = THETA*(1+0.61*mr_bfr);
    return THETAV

#==============================================================================

def mixratio2rh(pres,tair,mrat):    

   """ 
   This function estimates RH in % 
   From
   pres - Pressure in hPa
   tair - Temperature in °C
   mrat - Mixing ration in g/kg
   """
   eps       = 0.622 
   MRAT      = mrat/1000 
   vpres     = (MRAT*pres)/(eps+MRAT)
   sat_vpres = sat_vap_pres(tair)
   RH        = (vpres/sat_vpres)*100
   return RH    
   
#==============================================================================

def dewpt2rh(dpt,tair):
    """
    This function calculates Relative humidity (RH in %) 
    from Dewpoint temperature (dpt in °C) and 
    Air temperature (tair in °C)
    
    By Denny P. Alappattu
    Equation from http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
    """
    RH=100*(np.exp((17.625*dpt)/(243.04+dpt))/np.exp((17.625*tair)/(243.04+tair)))
    return RH
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   