#Qestion 8 : Le temps de calcul nécessaire pour simuler 10^6 simulation avec chacune des deux méthodes
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.special as ss
from matplotlib import rc
import time
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
nbs=1000000 #nombre de simulations
a=2.5 #alpha dans l'énoncé
lbCF=(a-1/(6*a))/(a-1) #lambda_CF dans l'énoncé

def estDansA(X) :    #Fonction qui teste si X appartient à A
    if ( (2/(a-1))*np.log(X[0])-np.log(lbCF*X[1]/X[0])+lbCF*X[1]/X[0]-1<=0 ) :
        return 1
    return 0

def estDansB(X) :    #Fonction qui teste si X appartient à B
    if ( (2/(a-1))*(X[0]-1)+X[0]/(lbCF*X[1])+lbCF*X[1]/X[0]-2<=0 ) :
        return 1
    return 0

def X_T1():  # La méthode du rejet usuelle
    X=[random(),random()]
    while(estDansA(X)!=1):
        X=[random(),random()]
    return X
def X_T2():  # La méthode du rejet en testant d'abord si (X_i £ B) puis si (X_i £ A)
    X=[random(),random()]
    while(estDansB(X)!=1 and estDansA(X)!=1):
        X=[random(),random()]
    return X
start_time=time.time()
Y=[X_T1() for i in range(nbs)]
print(time.time()-start_time)
start_time=time.time()
Y=[X_T2() for i in range(nbs)]
print(time.time()-start_time)