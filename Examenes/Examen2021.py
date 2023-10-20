# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:14:28 2023

@author: isidr
"""

from time import perf_counter
from matplotlib.pyplot import *
from numpy import *


mu = 2.0

def f(t, y):
    f1 = y[1]
    f2 = -mu*(y[0]**2-1)*y[1] + y[0]
    return(array([f1,f2]))


def rk12Sis(a, b, fun, y0, h0, tol):
    hmin = (b-a)*1.e-5  # paso de malla minimo
    hmax = (b-a)/1.  # paso de malla maximo

    # coeficientes RK
    q = 2  # orden del metodo mas uno
    A = zeros([q, q])
    A[1, 0] = 1.

    B = zeros(q)
    B[0] = 1.

    BB = zeros(q)
    BB[0] = 1./2.
    BB[1] = 1./2.

    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i, :])

    # inicializacion de variables
    t = array([a])  # nodos
    y = y0.reshape(len(y0), 1)   # soluciones
    h = array([h0])  # pasos de malla
    K = zeros([len(y0), q])
    k = 0  # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(q):
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]
                          * dot(A[i, :], transpose(K)))

        incrlow = dot(B, transpose(K))  # metodo de orden 4
        incrhigh = dot(BB, transpose(K))  # metodo de orden 5

        error = linalg.norm(h[k]*(incrhigh-incrlow), inf)  # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./5)  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)


# Datos del problema
a = 0 # extremo inferior del intervalo
b = 20 # extremo superior del intervalo
y0 = array([1.7,0.3]) # condicion inicial
tol = 1.e-6 # tolerancia
h0 = 0.1 # paso de malla inicial

tini = perf_counter()
(t1, y1, h1) = rk12Sis(a, b, f, y0, h0, tol) 
tfin = perf_counter()

figure('Ejercicio 1a')
subplot(311)
title('Graxica x frente a t')
plot(t1,y1[0]) 
# legend(['(1)','(3)'])

subplot(312)
title('Trayectoria')
plot(y1[0],y1[1]) 
# legend(['(1)','(3)'])

subplot(313)
title('Paso de h frente a t')
plot(t1, h1) 
# legend(['(1)','(3)'])

subplots_adjust(hspace=0.8)

gcf().suptitle("Ejercicio 1a")
show()

# Ejercicio 2

# def trapecioSis(a, b, fun, y0, h0, tol):
#     hmin = (b-a)*1.e-5  # paso de malla minimo
#     hmax = (b-a)/1.  # paso de malla maximo

#     # inicializacion de variables
#     t = array([a])  # nodos
#     y = y0.reshape(len(y0), 1)   # soluciones
#     h = array([h0])  # pasos de malla
#     k = 0  # contador de iteraciones

#     while (t[k] < b):
#         h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        
#         # Paso simple del método del trapecio
#         k1 = fun(t[k], y[:, k])
#         k2 = fun(t[k] + h[k], y[:, k] + h[k] * k1)
#         y1 = y[:, k] + h[k]/2 * (k1 + k2)
        
#         # Paso doble del método del trapecio
#         k1 = fun(t[k], y[:, k])
#         k2 = fun(t[k] + h[k]/2, y[:, k] + h[k]/2 * k1)
#         k3 = fun(t[k] + h[k], y[:, k] + h[k] * (-k1 + 2*k2))
#         y2 = y[:, k] + h[k]/6 * (k1 + 4*k2 + k3)

#         # Cálculo de la estimación del error
#         delta = linalg.norm(y2 - y1, ord=inf)/15

#         # Ajuste del tamaño del paso de integración
#         if delta < tol:
#             y = column_stack((y, y2))
#             t = append(t, t[k] + h[k])
#             hnew = min(h[k]*2, b-t[-1])
#         else:
#             hnew = max(h[k]/2, tol*abs(h[k]/delta)**0.25)*h[k] # fórmula de ajuste de paso de malla
#             # Asegurar que el siguiente paso no es menor que el paso de malla mínimo
#         if hnew < hmin:
#             hnew = hmin
    
#         # Asegurar que el siguiente paso no es mayor que el paso de malla máximo
#         if hnew > hmax:
#             hnew = hmax
    
#         h = append(h, hnew)  # actualización del paso de malla
#         k += 1  # actualizar el contador de iteraciones
    
#         # Solución encontrada
#         if t[k] == b:
#             break

#     return t, y





# # Datos del problema
# a = 0 # extremo inferior del intervalo
# b = 20 # extremo superior del intervalo
# y0 = array([1.7,0.3]) # condicion inicial
# tol = 1.e-6 # tolerancia
# h0 = 0.1 # paso de malla inicial

# tini = perf_counter()
# (t1, y1, h1) = trapecioSis(a, b, f, y0, h0, tol) 
# tfin = perf_counter()

# figure('Ejercicio 1a')
# subplot(311)
# title('Graxica x frente a t')
# plot(t1,y1[0]) 
# # legend(['(1)','(3)'])

# subplot(312)
# title('Trayectoria')
# plot(y1[0],y1[1]) 
# # legend(['(1)','(3)'])

# subplot(313)
# title('Paso de h frente a t')
# plot(t1, h1) 
# # legend(['(1)','(3)'])

# subplots_adjust(hspace=0.8)

# gcf().suptitle("Ejercicio 1a")
# show()








# Ejercicio 3
# El objetivo de este apartado es, para ε = 0.1, comparar la solución xε de (1) con condición inicial
# (3) con 2 + uε, donde uε es la solución de la ecuación (2) con condición inicial (4). Para ello,
# reescriba ambas ecuaciones diferenciales en forma equivalente como sistemas de primer orden
# y aplique el método RK4 con paso h = 0.05 para resolver las ecuaciones en el intervalo [0, 20].
# Ponga en pantalla el tiempo de cálculo correspondiente a uno y otro problema. Compare en
# una gráfica las aproximaciones obtenidas de xε(t) con las de 2 + uε(t). ¿Para qué valores de t se
# verifica (5)?