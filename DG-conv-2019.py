#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 17:24:01 2018

@author: gabrielbmiranda
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import csv
import copy as cp

omega = [0.,1.]

delta1 = -1/2
delta2 =  1/2

tau_u = 10
tau_p = 0

X = sp.symbols('x')
fonte = 4*(sp.pi**2)*sp.cos(2*sp.pi*X)
u_exa = 2*sp.pi*sp.sin(2*sp.pi*X)
p_exa = sp.cos(2*sp.pi*X)

def sol_analitica_u(pts):
    h = (omega[1]-omega[0])/(pts-1)
    x = []
    y = []
    for i in range(pts):
        x.append(omega[0]+i*h)
        y.append( 2*np.pi*np.sin(2*np.pi*x[i]) )
    return x, y

def sol_analitica_p(pts):
    h = (omega[1]-omega[0])/(pts-1)
    x = []
    y = []
    for i in range(pts):
        x.append(omega[0]+i*h)
        y.append( np.cos(2*np.pi*x[i]) )
    return x, y

w = []
t = []
def set_grau(grau):
    global w
    w = []
    global t
    t = []
    arquivo_w = open('weights.csv', newline='')
    leitor_w = csv.reader(arquivo_w)
    arquivo_t = open('abscissa.csv', newline='')
    leitor_t = csv.reader(arquivo_t)
    
    for row in leitor_w:
        if row[0] == 'quadrature_weights[' + str(grau+1) + ']':
            for i in range(1, grau+2):
                w.append(float(row[i]))
            break;
    
    for row in leitor_t:
        if row[0] == 'legendre_roots[' + str(grau+1) + ']':
            for i in range(1, grau+2):
                t.append(float(row[i]))
            break;
################################################

def gaussian_int(f,a,b,grau):
    fzao = lambda t : f.subs(X, ( (b-a)*t+b+a)/2 ) #altera f para fzao
    return (b-a)*sum([w[i]*fzao(t[i]) for i in range(grau+1)]) /2
    
def base_lagrange(x_,xi):
    base = 1
    for k in range(len(x_)):
        if x_[k] != xi:
            base *= (X-x_[k])/(xi-x_[k])
    return base 

def fim(v):
    return v.subs(X, omega[0] + h )

def ini(v):
    return v.subs(X, omega[0])

def matriz_local_A(grau,a,b):
    n = 2*(grau + 1)
    K = np.zeros((n,n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            K[2*i  ,2*j  ] = float( gaussian_int(  phi[j]* phi[i] + delta1* phi[j]* phi[i] + delta2*dphi[j]*dphi[i] ,a,b,grau) )
            K[2*i  ,2*j+1] = float( gaussian_int(- phi[j]*dphi[i] + delta1*dphi[j]* phi[i]                          ,a,b,grau) )
            K[2*i+1,2*j  ] = float( gaussian_int(-dphi[j]* phi[i] + delta1* phi[j]*dphi[i]                          ,a,b,grau) )
            K[2*i+1,2*j+1] = float( gaussian_int(                   delta1*dphi[j]*dphi[i]                          ,a,b,grau) )
    return K

def matriz_local_B(grau,a,b):
    n = 2*(grau + 1)
    K = np.zeros((n,n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            K[2*i  , 2*j ] = float( fim(phi[j]*phi[i])*tau_u )
            K[2*i+1,2*j+1] = float(-fim(phi[j]*phi[i])*tau_p )
            K[2*i  ,2*j+1] = float( fim(phi[j]*phi[i])       )
#            K[2*i+1,2*j  ] = float(-fim(phi[j]*phi[i])/2     )
    return K

def matriz_local_C(grau,a,b):
    n = 2*(grau + 1)
    K = np.zeros((n,n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            K[2*i  , 2*j ] = float( ini(phi[j]*phi[i])*tau_u )
            K[2*i+1,2*j+1] = float(-ini(phi[j]*phi[i])*tau_p )
            K[2*i  ,2*j+1] = float(-ini(phi[j]*phi[i])       )
#            K[2*i+1,2*j  ] = float( ini(phi[j]*phi[i])/2     )
    return K


def matriz_local_D(grau,a,b):
    n = 2*(grau + 1)
    K = np.zeros((n,n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            K[2*i  , 2*j ] = float(-ini(phi[j])*fim(phi[i])*tau_u )
            K[2*i+1,2*j+1] = float( ini(phi[j])*fim(phi[i])*tau_p )
#            K[2*i  ,2*j+1] = float(-ini(phi[j])*fim(phi[i])/2     )
#            K[2*i+1,2*j  ] = float(-ini(phi[j])*fim(phi[i])/2     )
    return K

def matriz_local_E(grau,a,b):
    n = 2*(grau + 1)
    K = np.zeros((n,n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            K[2*i  , 2*j ] = float(-fim(phi[j])*ini(phi[i])*tau_u )
            K[2*i+1,2*j+1] = float( fim(phi[j])*ini(phi[i])*tau_p )
#            K[2*i  ,2*j+1] = float(-fim(phi[j])*ini(phi[i])/2     )
#            K[2*i+1,2*j  ] = float( fim(phi[j])*ini(phi[i])/2     )
    return K


def vetor_fonte(grau,a,b,elemento):
    n = 2*(grau + 1)
    f = np.zeros(n)
    for j in range(int(n/2)):
        f[2*j  ] = float( gaussian_int(delta2*fonte*(dphi[j].subs(X, X - elemento*h)), a, b, grau) )
        f[2*j+1] = float( gaussian_int(-fonte*(phi[j].subs(X, X - elemento*h)),a,b,grau) )
    return f


def matriz_global(A,B,C,D,E,F0,FN,grau,nel):
    n = nel*(2*(grau+1))
    K = np.zeros((n,n))
    #Primeiro elemento da diagonal
    for i in range(2*(grau+1)):
        for j in range(2*(grau+1)):
            K[i,j] += A[i,j] + F0[i,j] + C[i,j]
    
    M = A + B + C
    for l in range(1,nel-1):
        for i in range(2*(grau+1)):
            for j in range(2*(grau+1)):
                K[l*(2*(grau+1))+i,l*(2*(grau+1))+j] += M[i,j]
    #Ultimo elemento da diagonal
    for i in range(2*(grau+1)):
        for j in range(2*(grau+1)):
            K[(nel-1)*(2*(grau+1))+i,(nel-1)*(2*(grau+1))+j] += A[i,j] + FN[i,j] + B[i,j]
    
    for l in range(nel-1):
        for i in range(2*(grau+1)):
            for j in range(2*(grau+1)):
                K[ (l+1)*2*(grau+1)+i,  l*2*(grau+1)+j ] += E[i,j]
                K[ l*2*(grau+1)+i,  (l+1)*2*(grau+1)+j ] += D[i,j]
    
    return K
    
def fonte_global(f,grau,nel):
    n = nel*(2*(grau+1))
    K = np.zeros(n)
    for i in range(nel):
        for j in range(2*(grau+1)):
            K[i*(2*(grau+1))+j] += f[i][j]
    return K

def contorno(grau):
    n = (2*(grau+1))
    K1 = np.zeros((n,n))
    K2 = np.zeros((n,n))
    for i in range(grau+1):
        for j in range(grau+1):
            K1[2*i  , 2*j ] = float( fim(phi[j]*phi[i])*tau_u )
            K1[2*i+1,2*j+1] = float(-fim(phi[j]*phi[i])*tau_p )
            K1[2*i  ,2*j+1] = float( fim(phi[j]*phi[i])       )
            K1[2*i+1,2*j+1] = float( 0 )
            
            K2[2*i  , 2*j ] = float( ini(phi[j]*phi[i])*tau_u )
            K2[2*i+1,2*j+1] = float(-ini(phi[j]*phi[i])*tau_p )
            K2[2*i  ,2*j+1] = float(-ini(phi[j]*phi[i])       )
            K2[2*i+1,2*j+1] = float( 0 )
    return K1, K2

def condicoes_contorno(K,F):
    F[-1] = 1
    F[ 1] = 1
    F[0] -= K[0,1]*F[1]
    for i in range( 2*(2*(grau+1)) ):
        F[i+2] -= K[i+2,1]*F[1]
        K[i+1,1] = 0.
        K[1,i+1] = 0.
        F[-(i+2)] -= K[-(i+2),-1]*F[-1]
        K[-1,-(i+2)] = 0.
        K[-(i+2),-1] = 0.
    K[-1,-1] = 1.
    K[ 1, 0] = 0.
    K[ 0, 1] = 0.
    K[ 1, 1] = 1.
    return K, F

pts, exata_u = sol_analitica_u(256)
pts, exata_p = sol_analitica_p(256)


refinamentos = 4

for k in range(1,2):
    grau = k
    
    erro_u = []
    erro_p = []
    tam_h = []
    for l in range(refinamentos):
        set_grau(k)
        nel = 2*2**(l+1)
        h = (omega[1]-omega[0])/nel
        
        tau_p = 0.000000000001* h
        tau_u = 1/tau_p
        
        tam_h.append(h)
        phi = [ ]        
        dphi = [ ]
        x_ = [ ]
        for i in range(grau+1):
            x_.append( omega[0]+i*h/grau )
        for i in range(grau+1):
            phi.append( base_lagrange(x_, x_[i]) )
            dphi.append( sp.diff(phi[-1],X))

        N = []
        A = matriz_local_A(grau, omega[0] ,omega[0]+h )
        B = matriz_local_B(grau, omega[0] ,omega[0]+h )
        C = matriz_local_C(grau, omega[0] ,omega[0]+h )
        D = matriz_local_D(grau, omega[0] ,omega[0]+h )
        E = matriz_local_E(grau, omega[0] ,omega[0]+h )
        
#        Condicoes de contorno
        F0, FN = contorno(grau)
        
        
        
        for i in range(nel):
            F = vetor_fonte(grau,omega[0]+i*h,omega[0]+(i+1)*h,i)
            N.append(F)
        
        a = matriz_global(A,B,C,D,E,F0,FN,grau,nel)
        b = fonte_global(N,grau,nel)
        
        a,b = condicoes_contorno(a,b)
        
        coeficientes = np.linalg.solve(a,b)
        
        
        u = []
        for j in range(nel*(grau+1)):
            u.append( coeficientes[2*j] ) 
        
        p = []
        for j in range(nel*(grau+1)):
            p.append( coeficientes[2*j+1] )
        
        locpts = []
        for i in range(grau+1):
            locpts.append(i*h/grau)
        
        temp = cp.copy(locpts)
        
        set_grau(grau+2)
        
        erroaux = float(gaussian_int( (sum( u[j]*phi[j] for j in range(grau+1) ) - u_exa )**2, locpts[0], locpts[-1], grau+2) )
        for i in range(nel-1):
            for j in range(grau+1):
                locpts[j] += h
            erroaux +=  float(gaussian_int( (sum( u[(i+1)*(grau+1)+j]*phi[j].subs(X, X - (i+1)*h) for j in range(grau+1) ) - u_exa )**2, locpts[0], locpts[-1], grau+2) )
        erro_u.append(np.sqrt(erroaux))
        
        locpts = temp
        erroaux = float(gaussian_int( (sum( p[j]*phi[j] for j in range(grau+1) ) - p_exa )**2, locpts[0], locpts[-1], grau+2) )
        for i in range(nel-1):
            for j in range(grau+1):
                locpts[j] += h
            erroaux +=  float(gaussian_int( (sum( p[(i+1)*(grau+1)+j]*phi[j].subs(X, X - (i+1)*h) for j in range(grau+1) ) - p_exa )**2, locpts[0], locpts[-1], grau+2) )
        erro_p.append(np.sqrt(erroaux))

    print('----------GRAU = ' + str(k) + '-----------')   
    print('----------pressao p----------')
#    for o in range(1, refinamentos):
#        print('Nel = ' + str(2*2**(o)) + ' p/ ' + str(3*2**(o+1)) +' & ' + str( (np.log(erro_p[o]) - np.log(erro_p[o-1]) ) / (np.log(tam_h[o])-np.log(tam_h[o-1])) ) + ' \\\ \hline')
    print(str( (np.log(erro_p[-1]) - np.log(erro_p[0]) ) / (np.log(tam_h[-1])-np.log(tam_h[0])) ))
    
    plt.figure(0)
    plt.plot(-np.log(tam_h), np.log(erro_p), '-x', label='k,l = ' + str(grau) )
    
    print('--------velocidade u---------')
#    for o in range(1, refinamentos):
#        print('Nel = ' + str(2*2**(o)) + ' p/ ' + str(2*2**(o+1)) +' & ' + str( (np.log(erro_u[o]) - np.log(erro_u[o-1]) ) / (np.log(tam_h[o])-np.log(tam_h[o-1])) ) + ' \\\ \hline')
    print(str( (np.log(erro_u[-1]) - np.log(erro_u[0]) ) / (np.log(tam_h[-1])-np.log(tam_h[0])) ))

    plt.figure(1)
    plt.plot(-np.log(tam_h), np.log(erro_u), '-x', label='k,l = ' + str(grau) )
    
#plt.title('Galerkin Descontínuo Misto - Taxa de Convergência p')
plt.grid()
plt.legend(loc='best')
plt.savefig('tax_conv-p_DG.eps')

plt.figure(0)
#plt.title('Galerkin Descontínuo Misto - Taxa de Convergência u')
plt.grid()
plt.legend(loc='best')
plt.savefig('tax_conv-u_DG.eps')