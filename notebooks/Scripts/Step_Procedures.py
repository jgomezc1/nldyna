# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 09:06:09 2019

@author: JULIAN PARRA
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def AnalyticalSLN(M,K,Xo,Vo,Po,W,DeltaT,Tt):
    """
    Funcion para determinar la respuesta dinamica analitica de un sistema de 1GDL sometido a
    vibración libre no amoritguada, ante una excitación armonica de la forma = PoSenWt
    
    INPUT:
    -----
    M  = Masa del sistema
    K  = Rigidez del sistema
    Xo = Depslazmaiento inicial
    Vo = Velocidad inicial
    Po = Amplitud de la excitación
    W  = Frecuencia de excitación 
    Tt = Tiempo total
    
    OUTPUT:
    -------
    U = Respuesta de desplazamientos analitica del sistema
    
    """
    
    Wn = np.sqrt(K/M)
    C  = (Po/K)*(1/(1-(W/Wn)**2))
    
    Ninc = int(round(Tt/DeltaT,0))
    U    = np.zeros((2,Ninc))
    
    for i in range (Ninc):
        t      = i*DeltaT 
        U[0,i] = t
        U[1,i] = C*(np.sin(W*t)-(W/Wn)*np.sin(Wn*t))
        
    #End for i    
         
    return U

def Force(Po,W,DeltaT,Tt):
    """
    INPUT:
    -----
    M  = Masa del sistema
    K  = Rigidez del sistema
    Xo = Depslazmaiento inicial
    Vo = Velocidad inicial
    Po = Amplitud de la excitación
    W  = Frecuencia de excitación 
    Tt = Tiempo total
    
    OUTPUT:
    -------
    F = Excitación del sistema
    
    """
    Ninc = int(round(Tt/DeltaT,0))
    P    = np.zeros((2,Ninc))
    
    for i in range (Ninc):
        t      = i*DeltaT 
        P[0,i] = t
        P[1,i] = Po*np.sin(W*t)
        
    #End for i   
    
    return P
    

def Increment_ThetaWilson(theta,dt,M,K,C,dtaF,U,V,A):
        
    c1 = 6.0/(theta*dt)
    c2 = theta*dt/2.0
    c3 = 3.0/(theta*dt)
    c4 = 6.0/(theta*theta*dt*dt)
    c5 = 1.0/theta
    c6 = dt/2.0
    c7 = dt*dt/2.0
    c8 = dt*dt/6.0
    
    a = c1*M + 3*C
    b = 3*M + c2*C
    
    dP   = theta*dtaF + a*V + b*A
    KE   = K + c3*C + c4*M
    dU   = dP/KE

    dA   = c4*dU - c1*V - 3*A
    dtaA = dA*c5
    dtaV = dt*A + dtaA*c6
    dtaU = dt*V + c7*A + c8*dtaA
    
    U = U + dtaU
    V = V + dtaV
    A = A + dtaA
    
    return U,V,A
 
def ThetaWilson_Resp(theta,DeltaT,M,K,C,Loading,Tt):
     
     Vo = 0
     Ao = 0
     Ninc = int(round(Tt/DeltaT,0))
     U = np.zeros(Ninc)
     T = np.zeros(Ninc)
     T[-1] = Ninc*DeltaT
     F = Loading[1,:]
     
     for k in range (Ninc-1):
          T[k] = k*DeltaT
          dtaF = F[k+1] - F[k] 
     
          U1, V1, A1 = Increment_ThetaWilson(theta,DeltaT,M,K,C,dtaF,U[k],Vo,Ao)
          
          U[k+1] = U1
          Vo = V1
          Ao = A1
     # End for k
     return U,T
  
def Plots(P, U1, T, U2):
     
     Fig = plt.figure(figsize=(15,4.5))
     Fig1 = Fig.add_subplot(121)
     Fig2 = Fig.add_subplot(122)
     
     Fig1.plot(P[0,:], P[1,:], '-.k')
     Fig1.set_title('Excitation')
     Fig1.set_xlabel('Time (seg)')
     Fig1.set_ylabel('Acceleration')
     Fig1.grid('On')
     
     
     Fig2.set_title('Displacement Response')
     Fig2.set_xlabel('Time (seg)')
     Fig2.set_ylabel('Displacement')
     AS = Fig2.plot(U2[0,:], U2[1,:],'gray',label='Analytical solution')
     TW = Fig2.plot(T, U1,'-.r', lw=2, label='The Theta Wilson method')
     plt.legend(bbox_to_anchor=(0.3, -0.7, -0.15, 0.5), borderaxespad=0.)
     Fig2.grid('on')      
     
     return Fig



    
    
    
    