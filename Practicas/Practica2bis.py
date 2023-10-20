# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 14:17:34 2023

@author: isidr
"""

from pylab import *
from time import perf_counter

# Para comprobar que las implementaciones de los métodos RK2(3), RK4(5) y RK5(4) son correctas, ejecuta los códigos para resolver el problema
# y = −y + 2 sin(t), t ∈ [0, 10],
# y(0) = π
# cuya solución exacta es la función y(t) = (π + 1)e^(-t) + sin(t) − cos(t). Los errores, número de iteraciones y pasos de malla 
# que se deben obtener (aproximadamente) son los siguientes:

# RK2(3): error = 4.4965 · 10^-5, iteraciones = 557, h mínimo = 0.005, h máximo = 0.1761.
# RK4(5): error = 7.4645 · 10^-6, iteraciones = 46, h mínimo = 0.005, h máximo = 0.4305.
# Se recomienda tomar h0 = (b − a)/2000, hmín = (b − a)/105, hmáx = (b − a)/10 y una tolerancia
#  para la estimación del error de tol = 10^-6. El error se calcula como el máximo entre los valores
#  absolutos de la diferencia entre la solución exacta y la solución aproximada en cada punto. 
#  Al calcular hmín y hmáx, no se deben tener en cuenta el último valor calculado (que no se utiliza) 
#  ni el penúltimo (que se ajusta para llegar exactamente al tiempo final). En el número de iteraciones 
#  no se debe tener en cuenta la condición inicial (es decir, es la longitud del array y menos 1).


def f(t, y):
    return -y + 2 * sin(t)

def exacta(t):
    return (pi + 1) * exp(-t) + sin(t) - cos(t)

def rk23(a, b, fun, y0, h0, tol):
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""
    # tol es la tolerancia de error

    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 3 # numero de etapas
    p = 2 # orden del método menos preciso
    A = zeros([q, q]) #coeficientes aij inicializado a 0 y ahora rellenamos
    A[1, 0] = 1./2.
    A[2, 0] = -1.
    A[2, 1] = 2. #tercera fila segunda columna
    
    B = zeros(q) #pesos
    B[1] = 1. #segundo (b2)
    
    BB = zeros(q) 
    BB[0] = 1./6.
    BB[1] = 2./3.
    BB[2] = 1./6.
    
    C = zeros(q) #No los pide, los calculas sumando las filas
    for i in range(q):
        C[i] = sum(A[i,:])
    # En un principio son arrays de longitud uno
    # inicializacion de variables
    t = array([a]) # nodos
    y = array([y0]) # soluciones
    h = array([h0]) # pasos de malla
    K = zeros(3) #array que contiene k1, k2, k3 que usa Runge-Kutta para avanzar
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q): #calcular K1, K2, K3
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i,:]*K))
        incrlow = sum(B*K) # metodo de orden 2 #phi
        # B es b1, b2, b3, K es k1, k2, k3 multiplica elemento a elemento y los suma
        incrhigh = sum(BB*K) # metodo de orden 3 #phiEstrella
        # BB es b1*, b2*, b3*, K es k1, k2, k3 multiplica elemento a elemento y los suma
        error = h[k]*(incrhigh-incrlow) # estimacion del error
        y = append(y, y[k]+h[k]*incrlow) # y_(k+1)
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1)) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    
    print("Numero de iteraciones: " + str(k-1))
    return (t, y, h)


print('EJERCICIO 1a')
    
# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = pi # condicion inicial
h0 = (b-a)/2000. #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rk23(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
figure('Ej1')
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK2(3)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')


def rk45(a, b, fun, y0, h0, tol): #de quinto orden necesita al menos 6 etapas
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""
    # tol es la tolerancia de error

    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 6 # numero de etapas
    p = 4 # orden del método menos preciso
    A = zeros([q, q]) #coeficientes aij inicializado a 0 y ahora rellenamos
    A[1, 0] = 1./4.
    A[2, 0] = 3./32.
    A[2, 1] = 9./32.
    A[3, 0] = 1932./2197.
    A[3, 1] = -7200./2197.
    A[3, 2] = 7296./2197.
    A[4, 0] = 439./216.
    A[4, 1] = -8.
    A[4, 2] = 3680./513.
    A[4, 3] = -845./4104.
    A[5, 0] = -8./27.
    A[5, 1] = 2.
    A[5, 2] = -3544./2565.
    A[5, 3] = 1859./4104.
    A[5, 4] = -11./40.

    
    B = zeros(q) #pesos
    B[0] = 25./216.
    B[2] = 1408./2565.
    B[3] = 2197./4104.
    B[4] = -1./5.
    
    BB = zeros(q) 
    BB[0] = 16./135.
    BB[2] = 6656./12825.
    BB[3] = 28561./56430.
    BB[4] = -9./50.
    BB[5] = 2./55.

    
    C = zeros(q) #No los pide, los calculas sumando las filas
    for i in range(q):
        C[i] = sum(A[i,:])
    # En un principio son arrays de longitud uno
    # inicializacion de variables
    t = array([a]) # nodos
    y = array([y0]) # soluciones
    h = array([h0]) # pasos de malla
    K = zeros(6) #array que contiene k1, k2, k3 que usa Runge-Kutta para avanzar
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q): #calcular K1, K2, K3
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i,:]*K))
        incrlow = sum(B*K) # metodo de orden 2 #phi
        # B es b1, b2, b3, K es k1, k2, k3 multiplica elemento a elemento y los suma
        incrhigh = sum(BB*K) # metodo de orden 3 #phiEstrella
        # BB es b1*, b2*, b3*, K es k1, k2, k3 multiplica elemento a elemento y los suma
        error = h[k]*(incrhigh-incrlow) # estimacion del error
        y = append(y, y[k]+h[k]*incrlow) # y_(k+1)
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1)) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    
    print("Numero de iteraciones: " + str(k-1))
    return (t, y, h)

print('EJERCICIO 1b')

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 30. # extremo superior del intervalo
y0 = 0. # condicion inicial
h0 = 0.1 #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rk45(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
figure('Ej2')
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK4(5)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')











def f(t,y):
    f1= 3*y[0] - 2*y[1]
    f2= - y[0] + 3*y[1] - 2*y[2]
    f3= - 2*y[1] + 3*y[2]
    return array([f1,f2, f3])

def exacta(t):
    f1= -exp(5*t)/4. + 3*exp(3*t)/2. - exp(t)/4.
    f2= exp(5*t)/4. - exp(t)/4.
    f3= -exp(5*t)/8. - 3*exp(3*t)/4. - exp(t)/8.
    return(array([f1,f2, f3]))


def rkSistemas23(a, b, fun, y0, h0, tol):
    
    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 3 # orden del metodo mas uno
    A = zeros([q, q])
    A[1, 0] = 1./2.
    A[2, 0] = -1.
    A[2, 1] = 2.
    
    B = zeros(q)
    B[1] = 1.
    
    BB = zeros(q)
    BB[0] = 1./6.
    BB[1] = 2./3.
    BB[2] = 1./6.
    
    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i,:])
    
    # inicializacion de variables
    t = array([a]) # nodos
    # y = y0 # soluciones
    y = y0.reshape(len(y0), 1)
    h = array([h0]) # pasos de malla
    K = zeros([len(y0),q])
    k = 0 # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(3):
            K[:,i] = fun(t[k]+C[i]*h[k], y[:,k]+h[k]*dot(A[i,:],transpose(K)))
        
        incrlow = dot(B,transpose(K)) # metodo de orden 2
        incrhigh = dot(BB,transpose(K)) # metodo de orden 3
            
        error = norm(h[k]*(incrhigh-incrlow),inf) # estimacion del error
        y = column_stack((y, y[:,k]+h[k]*incrlow))
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./q) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    print("Numero de iteraciones: " + str(k-1))   
    return (t, y, h)

print('EJERCICIO 2a')
    
# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = array([1, 0, -1]) # condicion inicial
h0 = (b-a)/2000. #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rkSistemas23(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
figure('Ej2a')
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK2(3)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')


def rkSistemas45(a, b, fun, y0, h0, tol):
    
    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 6 # orden del metodo mas uno
    A = zeros([q, q])
    A[1, 0] = 1./4.
    A[2, 0] = 3./32.
    A[2, 1] = 9./32.
    A[3, 0] = 1932./2197.
    A[3, 1] = -7200./2197.
    A[3, 2] = 7296./2197.
    A[4, 0] = 439./216.
    A[4, 1] = -8.
    A[4, 2] = 3680./513.
    A[4, 3] = -845./4104.
    A[5, 0] = -8./27.
    A[5, 1] = 2.
    A[5, 2] = -3544./2565.
    A[5, 3] = 1859./4104.
    A[5, 4] = -11./40.
    
    B = zeros(q)
    B[0] = 25./216.
    B[2] = 1408./2565.
    B[3] = 2197./4104.
    B[4] = -1./5.
    
    BB = zeros(q)
    BB[0] = 16./135.
    BB[2] = 6656./12825.
    BB[3] = 28561./56430.
    BB[4] = -9./50.
    BB[5] = 2./55.
    
    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i,:])
    
    # inicializacion de variables
    t = array([a]) # nodos
    # y = y0 # soluciones
    y = y0.reshape(len(y0), 1)
    h = array([h0]) # pasos de malla
    K = zeros([len(y0),q])
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q):
            K[:,i] = fun(t[k]+C[i]*h[k], y[:,k]+h[k]*dot(A[i,:],transpose(K)))
        
        incrlow = dot(B,transpose(K)) # metodo de orden 4
        incrhigh = dot(BB,transpose(K)) # metodo de orden 5
            
        error = norm(h[k]*(incrhigh-incrlow),inf) # estimacion del error
        y = column_stack((y, y[:,k]+h[k]*incrlow))
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./5) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    print("Numero de iteraciones: " + str(k-1))        
    return (t, y, h)

print('EJERCICIO 2b')

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = array([1, 0, -1]) # condicion inicial
h0 = (b-a)/2000. #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rkSistemas45(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
figure('Ej2b')
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK4(5)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')














