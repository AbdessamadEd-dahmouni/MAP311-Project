#Question 7 : Les probabilités de rejet par la méthode de Monte-Carlo 
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
C=np.exp(a-1)*ss.gamma(a)/(2*np.power((a-1),a))
pRejetCF=1-C/lbCF   #la valeur théorique de la probabilité de rejet pour lamba_CF
pRejet_s=1-C/lb_s   #la valeur héorique de la probabilité de rejet pour lambda soulignée

def estDansA(X,lb) :    #teste si X appartient à A
    if ( (2/(a-1))*np.log(X[0])-np.log(lb*X[1]/X[0])+lb*X[1]/X[0]-1<=0 ) :
        return 1
    return 0
    
N=0
M=0
X=[random(),random()]
for i in range(nbs):
    if estDansA(X,lb_s)==0 :
        N+=1
    if estDansA(X,lbCF)==0 :
        M+=1
    X=[random(),random()]
print("la probabilité de rejet pour _lambda_ par la méthode de Monte-Carlo est : "+str(N/nbs)+", la valeur théorique est : "+str(pRejet_s))
print("la probabilité de rejet pour lambdaCF par la méthode de Monte-Carlo est :"+str(M/nbs)+", la valeur théorique est : "+str(pRejetCF))