# -*- coding: utf-8 -*-
"""
Created on Sun May 14 09:27:45 2023

@author: isidr
"""

from numpy import *
from matplotlib.pyplot import *
from time import perf_counter
    

def AB3(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])

    for k in range(2):
        k1 = f[k]
        k2 = fun(t[k]+0.5*h, y[k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[k]+h*k3)
        y[k+1] = y[k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[k+1] = fun(t[k+1],y[k+1])
    
    for k in range(2, N):
        y[k+1] = y[k] + h/12 * (23*f[k] - 16*f[k-1] + 5*f[k-2])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)


# Ejercicio 1


def fun(t, y):
    return (-y + 2*sin(t))


def exacta(t):
    return ((pi + 1)*exp(-t) + sin(t) - cos(t))


y0 = pi
a = 0.0
b = 10

N = 50


tini = perf_counter()
(t, y) = AB3(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('AB3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

def AM3_punto_fijo(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[k]
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        z = y[k]
        contador = 0
        distancia = 1 + 1e-12
        C = y[k] + h/24 * (19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = 9*h*fun(t[k+1], z)/24 + C
            distancia = abs(z - z_nuevo)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[k+1] = z
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y, iteraciones)

tini = perf_counter()
(t, y, maxiter) = AM3_punto_fijo(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('AM3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('maxiter : '+str(maxiter))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

def ABM3(a,b,fun, N, y0):
    
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a)/float(N) 
    y[0] = y0
    f[0] = fun(a,y[0])

    for k in range(2):  ### Arranco con RK4 estandar
        k1 = f[k]
        k2 = fun(t[k]+0.5*h, y[k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[k]+h*k3)
        y[k+1] = y[k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[k+1] = fun(t[k+1],y[k+1])
    for k in range(2,N):
        Ck = y[k] +h/24*(19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k]+h
       
        
        z0 = y[k]+h/12*(23*f[k] - 16*f[k-1] + 5*f[k-2])
        

        y[k+1] = h*9/24*fun(t[k+1],z0) + Ck
        f[k+1] = fun(t[k+1],y[k+1])
    

    return(t,y)

tini = perf_counter()
(t, y) = ABM3(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('ABM3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

legend(('AB3', 'AM3', 'ABM3'))
show()


# Ejercicio 2


def fun2(t, y):
    f1 = 3*y[0] - 2*y[1]
    f2 = -y[0] + 3*y[1] - 2*y[2]
    f3 = -y[1] + 3*y[2]
    return array([f1, f2, f3])


def exacta2(t):
    f1 = -1/4 * exp(5*t) + 3/2 * exp(3*t) - 1/4 * exp(t)
    f2 = 1/4 * exp(5*t) - 1/4 * exp(t)
    f3 = -1/8 * exp(5*t) - 3/4 * exp(3*t) - 1/8 * exp(t)
    return (array([f1, f2, f3]))


y0 = array([1, 0, -1])

a = 0.0
b = 1

N = 100

def AB3_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])

    for k in range(2):
        y[:, k+1] = y[:, k] + h*fun(t[k] + h/2, y[:, k] + h/2 * f[:, k])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        y[:, k+1] = y[:, k] + h/12 * (23*f[:, k] - 16*f[:, k-1] + 5*f[:, k-2])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
        
    return (t, y)

tini = perf_counter()
(t, y) = AB3_sistema(a, b, fun2, N, y0)
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('AB3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

def AM3pfSis(a,b,fun,N,y0):
    
    tol = 1.e-12
    stop = 200
    
    y = zeros([len(y0),N+1])
    t = zeros(N+1)
    f = zeros([len(y0),N+1])
    t[0] = a
    h = (b-a)/float(N) 
    y[:,0] = y0
    f[:,0] = fun(a,y0)
    maxiter = 0
    for k in range(2):  ### Arranco con RK4 estandar
        k1 = fun(t[k],y[:,k])
        k2 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[:,k]+h*k3)
        y[:,k+1] = y[:,k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    for k in range(2,N):
        Ck = y[:,k] +h/24*(19*f[:,k] - 5*f[:,k-1] + f[:,k-2])
        t[k+1] = t[k]+h
        cont = 0
        z = y[:,k]
        error = tol +1
        while(error >=tol and cont < stop):
            znew = h*9/24*fun(t[k+1], z) +Ck
            error = max(abs(znew - z))
            z = znew
            cont += 1
        maxiter = max(maxiter,cont)
        if (cont == stop):
            print('El método de punto fijo no ha convergido tras '+str(stop)+' iteraciones')
        y[:,k+1] = z
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    return(t,y,maxiter)


tini = perf_counter()
(t, y, maxiter) = AM3pfSis(a, b, fun2, N, y0)
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('AM3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Maxiter = '+str(maxiter))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

def ABM3Sis(a,b, fun, N,y0):
    
    tol = 1.e-12
    stop = 200
    
    y = zeros([len(y0),N+1])
    t = zeros(N+1)
    f = zeros([len(y0),N+1])
    t[0] = a
    h = (b-a)/float(N) 
    y[:,0] = y0
    f[:,0] = fun(a,y0)

    for k in range(2):  ### Arranco con RK4 estandar
        k1 = fun(t[k],y[:,k])
        k2 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[:,k]+h*k3)
        y[:,k+1] = y[:,k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    for k in range(2,N):
        Ck = y[:,k] +h/24*(19*f[:,k] - 5*f[:,k-1] + f[:,k-2])
        t[k+1] = t[k]+h
       
        
        z0 = y[:,k]+h/12*(23*f[:,k] - 16*f[:,k-1] + 5*f[:,k-2])
        

        y[:,k+1] = h*9/24*fun(t[k+1],z0) + Ck
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    

    return(t,y)


tini = perf_counter()
(t, y) = ABM3Sis(a, b, fun2, N, y0)
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('ABM3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Tiempo CPU= '+str(tcpu))
print('---------------')