# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:14:28 2023

@author: isidr
"""

from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
# 1. Se considera el problema de Cauchy
# y' = 4 − 1000y,
# y(0) = 10.
# Dibuje la frontera de la región de estabilidad absoluta del método del punto medio, cuyo tablero de Butcher es
# 0      0       0
# 1/2   1/2      0
#        0       1
# Estime el menor número de particiones N del intervalo [0, 1] que hay que tomar para que no se
# produzcan oscilaciones no controladas al aplicar este método al problema (1). Compruebe que la
# estimación es correcta aplicando el método del punto medio con valores ligeramente mayores y
# menores que el N estimado.
def f(t,y):
    return 4 - 1000*y

def exacta(t):
    return  4/1000 + (10-4/1000)*exp(-1000*t)

# def rk(a, b, fun, N, y0):
#     """Implementacion del metodo del punto medio en el intervalo [a, b]"""
#     h = (b-a)/N  # paso de malla
#     t = zeros(N+1)  # inicializacion del vector de nodos
#     y = zeros(N+1)  # inicializacion del vector de resultados
#     t[0] = a  # nodo inicial
#     y[0] = y0  # valor inicial

#     # Metodo del punto medio
#     for k in range(N):
#         k1 = fun(t[k], y[k])
#         k2 = fun(t[k] + h/2, y[k] + h*k1/2)
#         t[k+1] = t[k] + h
#         y[k+1] = y[k] + h*k2
#     return (t, y)
# #Datos del problema
# a = 0
# b = 1
# y0 = 10
# malla = [400, 450, 500, 550, 600]

# for i in malla:
#     (t, y) = rk(a, b, f, i, y0)
#     ye = exacta(t)
#     error = abs(ye - y)
#     print('N ='+ str(i), 'error =', error[-1])
#     plot(t, y, '-*')
#     plot(t, ye, 'k')
#     xlabel('t')
#     ylabel('y')
#     legend(['RK N=' + str(i), 'exacta'])
#     grid(True)
#     show()


# Ejercicio 2
# Se considera el método implícito cuyo tablero de Butcher es
# 1/2   1/2
#        1 
# a) Esciba un program python que resuelva la ecuación diferencial (1) usando este método.
# Indicación. Observe que, en este caso, la ecuación a resolver para calcular la etapa y(1)
# k es lineal: despéjela y use las expresión hallada para programar el método.
# b) Aplique el programa realizado al problema (1) en el intervalo [0, 1] tomando particiones
# uniformes de N = 20, 40, 80, 160, 320 subintervalos. Calcule y ponga en pantalla los tiempos
# de cálculo y los errores cometidos. Compare la gráfica de la solución exacta con las distintas
# aproximaciones numéricas obtenidas.

def rk(a, b, fun, N, y0):
    """Implementacion del metodo del punto medio en el intervalo [a, b]"""
    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo del punto medio
    for k in range(N):
        k1 = fun(t[k] + h/2, (y[k] + 2*h)/(1 + 500*h))
        t[k+1] = t[k] + h*1/2
        y[k+1] = y[k] + h*k1
    return (t, y)

#Datos del problema
a = 0
b = 1
y0 = 10
malla = [20, 40, 80, 160, 320]

for i in malla:
    tin = perf_counter()
    (t, y) = rk(a, b, f, i, y0)
    tout = perf_counter()
    ye = exacta(t)
    error = abs(ye - y)
    print('N ='+ str(i), 'error =', error[-1])
    print('Tiempo =', tout - tin)
    plot(t, y, '-*')
    plot(t, ye, 'k')
    xlabel('t')
    ylabel('y')
    if i == 20:
        leyenda = ['RK N = 40']
    if i > 20:
        leyenda.append('RK N = '+ str(i))
    if i == 320:
        leyenda.append('Exacta')
        legend(leyenda)
        grid(True)
show()