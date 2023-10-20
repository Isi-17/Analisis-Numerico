# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *

def euler(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h*fun(t[k], y[k])
    
    return (t, y)

def taylor2(a, b, f, df, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0
    
    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h*f(t[k], y[k]) + h**2/2*df(t[k], y[k])
        
    return (t, y)

def taylor3(a, b, f, d1f, d2f, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0
    
    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h*f(t[k], y[k]) + h**2/2*d1f(t[k], y[k]) + h**3/6*d2f(t[k], y[k])
        
    return (t, y)

def heun(a, b, fun, N, y0):    
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h/2 * (fun(t[k], y[k]) + fun(t[k+1], y[k] + h*fun(t[k], y[k])))
    
    return (t, y)

def punto_medio(a, b, fun, N, y0):    
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h * fun(t[k] + h/2, y[k] + h/2 * fun(t[k], y[k]))
    
    return (t, y)

def rk4(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a 
    y[0] = y0 

    for k in range(N):
        t[k+1] = t[k] + h
        k1 = fun(t[k], y[k])
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
    
    return (t, y)

def euler_sistema(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0), N+1])
    t[0] = a
    y[:, 0] = y0

    for k in range(N):
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k]+h*fun(t[k], y[:, k])
    
    return (t, y)
    
def heun_sistema(a, b, fun, N, y0):    
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0), N+1])
    t[0] = a
    y[:, 0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[:, k+1] = y[:, k] + h/2 * (fun(t[k], y[:, k]) + fun(t[k+1], y[:, k] + h*fun(t[k], y[:, k])) )
    
    return (t, y)

def punto_medio_sistema(a, b, fun, N, y0):    
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0), N+1])
    t[0] = a
    y[:, 0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[:, k+1] = y[:, k] + h * fun(t[k] + h/2, y[:, k] + h/2 * fun(t[k], y[:, k]))
    
    return (t, y)

def rk4_sistema(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0), N+1])
    t[0] = a 
    y[:, 0] = y0

    for k in range(N):
        t[k+1] = t[k]+h
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k] + h/2, y[:, k] + h/2 *k1)
        k3 = fun(t[k] + h/2*h, y[:, k] + h/2 *k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
    
    return (t, y)

def rk23(a, b, fun, y0, h0, tolerancia):
    hmin, hmax = (b-a)*1e-5, (b-a)/10
    
    A = array([[0, 0, 0], [1/2, 0, 0], [-1, 2, 0]])
    B = array([0, 1, 0])
    BB = array([1/6, 2/3, 1/6])
    C = array([0, 1/2, 1])
    
    t = array([a])
    y = array([y0])
    h = array([h0])
    k = 0
    
    while (t[k] < b):
        h[k] = min(h[k], b - t[k])
        K = zeros(3)
        for i in range(3):
            K[i] = fun(t[k] + C[i]*h[k], y[k] + h[k]*sum(A[i,:]*K))
        incremento_2, incremento_3 = sum(B*K), sum(BB*K)
        error = h[k] * (incremento_3 - incremento_2)
        y = append(y, y[k] + h[k]*incremento_2)
        t = append(t, t[k] + h[k])
        h_nuevo = min(max(0.9 * h[k] * abs(tolerancia/error)**(1/3), hmin), hmax)
        h = append(h, h_nuevo)
        k += 1

    return (t, y, h)

def rk45(a, b, fun, y0, h0, tol):  
    hmin, hmax = (b-a)*1e-5, (b-a)/10
    
    q = 6
    A = zeros([q, q])
    A[1, 0] = 1/4
    A[2, 0] = 3/32
    A[2, 1] = 9/32
    A[3, 0] = 1932/2197
    A[3, 1] = -7200/2197
    A[3, 2] = 7296/2197
    A[4, 0] = 439/216
    A[4, 1] = -8
    A[4, 2] = 3680/513
    A[4, 3] = -845/4104
    A[5, 0] = -8/27
    A[5, 1] = 2
    A[5, 2] = -3544/2565
    A[5, 3] = 1859/4104
    A[5, 4] = -11/40
    
    B = zeros(q)
    B[0] = 25/216
    B[2] = 1408/2565
    B[3] = 2197/4104
    B[4] = -1/5
    
    BB = zeros(q)
    BB[0] = 16/135
    BB[2] = 6656/12825
    BB[3] = 28561/56430
    BB[4] = -9/50
    BB[5] = 2/55
    
    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i, :])
    
    t = array([a])
    y = array([y0])
    h = array([h0])
    k = 0
    
    while (t[k] < b):
        h[k] = min(h[k], b - t[k])
        K = zeros(q)
        for i in range(q):
            K[i] = fun(t[k] + C[i]*h[k], y[k] + h[k]*sum(A[i,:]*K))
        incrlow = sum(B*K)
        incrhigh = sum(BB*K)
        error = h[k] * (incrhigh - incrlow)
        y = append(y, y[k] + h[k]*incrlow)
        t = append(t, t[k] + h[k])
        hnew = 0.9 * h[k] * abs(tol/error)**(1/5)
        hnew = min(max(hnew, hmin), hmax)
        h = append(h, hnew)
        k += 1

    return (t, y, h)

def rk23_sistema(a, b, fun, y0, h0, tolerancia):
    hmin, hmax = (b-a)*1e-5, (b-a)/10
    
    A = array([[0, 0, 0], [1/2, 0, 0], [-1, 2, 0]])
    B = array([0, 1, 0])
    BB = array([1/6, 2/3, 1/6])
    C = array([0, 1/2, 1])
    
    t = array([a])
    y =zeros([len(y0), 1])
    y[:, 0] = y0
    h = array([h0])
    K = zeros([len(y0), 3])
    k = 0
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k])
        for i in range(3):
            K[:,i] = fun(t[k] + C[i]*h[k], y[:, k] + h[k]*dot(A[i, :], transpose(K)))
        incremento_2, incremento_3 = dot(B, transpose(K)), dot(BB, transpose(K))
        error = norm(h[k]*(incremento_3 - incremento_2), inf)
        y = column_stack((y, y[:, k] + h[k]*incremento_2))
        t = append(t, t[k] + h[k])
        h_nuevo = min(max(0.9 * h[k] * abs(tolerancia/error)**(1/3), hmin), hmax)
        h = append(h, h_nuevo)
        k += 1
    
    return (t, y, h)

def rk45_sistema(a, b, fun, y0, h0, tol):  
    hmin, hmax = (b-a)*1e-5, (b-a)/10
    
    q = 6
    A = zeros([q, q])
    A[1, 0] = 1/4
    A[2, 0] = 3/32
    A[2, 1] = 9/32
    A[3, 0] = 1932/2197
    A[3, 1] = -7200/2197
    A[3, 2] = 7296/2197
    A[4, 0] = 439/216
    A[4, 1] = -8
    A[4, 2] = 3680/513
    A[4, 3] = -845/4104
    A[5, 0] = -8/27
    A[5, 1] = 2
    A[5, 2] = -3544/2565
    A[5, 3] = 1859/4104
    A[5, 4] = -11/40
    
    B = zeros(q)
    B[0] = 25/216
    B[2] = 1408/2565
    B[3] = 2197/4104
    B[4] = -1/5
    
    BB = zeros(q)
    BB[0] = 16/135
    BB[2] = 6656/12825
    BB[3] = 28561/56430
    BB[4] = -9/50
    BB[5] = 2/55
    
    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i, :])
    
    t = array([a])
    y = zeros([len(y0), 1])
    y[:,0] = y0
    h = array([h0])
    K = zeros([len(y0), q])
    k = 0
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k])
        for i in range(q):
            K[:,i] = fun(t[k]+C[i]*h[k], y[:,k]+h[k]*dot(A[i,:],transpose(K)))
        incrlow = dot(B, transpose(K))
        incrhigh = dot(BB, transpose(K))
        error = norm(h[k]*(incrhigh-incrlow),inf)
        y = column_stack((y, y[:,k] + h[k]*incrlow))
        t = append(t, t[k]+h[k])
        hnew = 0.9*h[k]*abs(tol/error)**(1/5)
        hnew = min(max(hnew,hmin),hmax)
        h = append(h, hnew)
        k += 1

    return (t, y, h)


def AB2(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N 
    y[0] = y0
    f[0] = fun(a, y[0])
    y[1] = y[0] + h*f[0]
    t[1] = a + h
    f[1] = fun(t[1], y[1])

    for k in range(1, N):
        y[k+1] = y[k] + h/2 * (3*f[k] - f[k-1])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)

def AB3(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])

    for k in range(2):
        y[k+1] = y[k] + h*fun(t[k] + h/2, y[k] + h/2 * f[k])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        y[k+1] = y[k] + h/12 * (23*f[k] - 16*f[k-1] + 5*f[k-2])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)

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

def AM3_newton(a, b, fun, dfun, N, y0):
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
            z_nuevo = z - (z - 9*h*fun(t[k+1], z)/24 - C)/(1 - 9*h*dfun(t[k+1], z)/24)
            distancia = abs(z - z_nuevo)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[k+1] = z
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y, iteraciones)

def predictor_corrector(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
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
        t[k+1] = t[k] + h
        y_nuevo = y[k] + h/12 * (23*f[k] - 16*f[k-1] + 5*f[k-2])
        f_nuevo = fun(t[k+1], y_nuevo)
        y[k+1] = y[k] + h/24 * (9*f_nuevo + 19*f[k] - 5*f[k-1] + f[k-2])
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y)

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


def AM3_punto_fijo_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[:, k]
        k2 = fun(t[k] + h/2, y[:, k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[:, k] + h/2 * k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        z = y[:, k]
        contador = 0
        distancia = 1 + 1e-12
        C = y[:, k] + h/24 * (19*f[:, k] - 5*f[:, k-1] + f[:, k-2])
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = 9*h*fun(t[k+1], z)/24 + C
            distancia = max(abs(z - z_nuevo))
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[:, k+1] = z
        f[:, k+1] = fun(t[k+1], y[:, k+1])

    return (t, y, iteraciones)

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

def predictor_corrector_sistema(a, b, fun, N, y0):  #ABM3Sis
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[:, k]
        k2 = fun(t[k] + h/2, y[:, k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[:, k] + h/2 * k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        t[k+1] = t[k] + h
        y_nuevo = y[:, k] + h/12 * (23*f[:, k] - 16*f[:, k-1] + 5*f[:, k-2])
        f_nuevo = fun(t[k+1], y_nuevo)
        y[:, k+1] = y[:, k] + h/24 * (9*f_nuevo + 19*f[:, k] - 5*f[:, k-1] + f[:, k-2])
        f[:, k+1] = fun(t[k+1], y[:, k+1])

    return (t, y)

def localizar_frontera_RK(dR, N):
    Npoints = 5000
    h = 2*N*pi/Npoints
    z = zeros(Npoints +1 , dtype = complex)
    z[0] = 0
    t = 0
    for k in range(len(z)-1):
        z[k+1] = z[k]+ h*1j*exp(1j*t)/dR(z[k])
        t = t + h
    x, y = real(z), imag(z)
    plot(x, y)
    grid(True)
    axis('equal')

def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')
