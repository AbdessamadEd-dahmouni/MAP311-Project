#Question 8 : Intervalle de confiance de la probabilité P(X£B)
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.special as ss
import sympy as sy
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
nbs=10000 #nombre de simulations
a=2.5 #alpha dans l'énoncé
lbCF=(a-1/(6*a))/(a-1) #lambda_CF dans l'énoncé
lb_s=np.power(((a+1)/(a-1)),(a+1)/2) #lambda soulignée dans l'énoncé

def estDansB(X,lb) :    #teste si X appartient à A
    if ( (2/(a-1))*(X[0]-1)+X[0]/(lb*X[1])+lb*X[1]/X[0]-2<=0 ) :
        return 1
    return 0

n=0
m=0
X=[random(),random()]
for i in range(nbs):
    if estDansB(X,lb_s) :
        n+=1
    if estDansB(X,lbCF) :
        m+=1
    X=[random(),random()]
    
print("l'intervalle de confiance pour _lambda_ est ["+str(n/nbs-1/np.sqrt(nbs))+","+str(n/nbs+1/np.sqrt(nbs))+"]")
print("l'intervalle de confiance pour lambdaCF est ["+str(m/nbs-1/np.sqrt(nbs))+","+str(m/nbs+1/np.sqrt(nbs))+"]")