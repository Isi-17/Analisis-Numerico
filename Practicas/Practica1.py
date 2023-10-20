# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 12:09:15 2020

@author: Usuario
"""

from pylab import *
from time import perf_counter
# Mide el tiempo CPU que tarda en metodos numericos

def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 0.5*(t**2 - y)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return t**2 - 4*t + 8 - 7.*exp(-0.5*t)

def euler(a, b, fun, N, y0):
    """Implementacion del metodo de Euler en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k])
        t[k+1] = t[k]+h
    
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
N = 20 # numero de particiones
y0 = 1. # condicion inicial

tini = perf_counter() # Tiempo CPU antes de llamar al metodo

(t, y) = euler(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter() # Tiempo CPU despues de llamar al metodo

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
figure('Ej0')
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Euler', 'exacta'])
grid(True)

# Calculo del error cometido
error = max(abs(y-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

show() # muestra la grafica

print('EJERCICIO 1a')
# Compruebe que, para el problema considerado, el orden del método de Euler es uno. Para
#ello, realice una serie de ejecuciones del programa con N = 10, 20, 40, 80, 160. Obtenga en
#cada caso el error eN que se comete y calcule los cocientes eN/e2N, N =10, 20, 40, 80. Dibuje
#en una gráfica las aproximaciones obtenidas con los distintos valores de N junto con la
#solución exacta

a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
N = 20 # numero de particiones
y0 = 1. # condicion inicial
figure('Ej1a')
for k in [10, 20, 40, 80, 160]:
    (t, y) = euler(a, b, f, k, y0) # llamada al metodo de Euler
    ye = exacta(t) # calculo de la solucion exacta
    # Dibujamos las soluciones
    plot(t, y, '-*') # dibuja la solucion aproximada
    if k == 10:
        leyenda = ['Euler N = 10']
    if k > 10:
        leyenda.append('Euler N = '+ str(k))
    
    # Calculo del error cometido
    error = max(abs(y-ye))
    # Resultados
    print('Error',k,': ' + str(error))
    if k!= 160:
        (ts, ys) = euler(a, b, f, 2*k, y0) # llamada al metodo de Euler
        yes = exacta(ts) # calculo de la solucion exacta
        errors = max(abs(ys-yes))
        print('Error e',k,'/e',2*k,': ' + str(error/errors))
    print('-----')
    if k == 160:
        leyenda.append('Exacta')
        plot(t, ye, 'k') # dibuja la solucion exacta
        xlabel('t')
        ylabel('y')
        legend(leyenda)
        grid(True)

show() # muestra la grafica
print('EJERCICIO 1b')
# En un depósito de 20 l, que inicialmente está lleno de agua pura, empieza a entrar, desde el
#instante t = 0, agua con una concentración de 3 grs/l de sal a razón de 2 l por segundo. Al
#mismo tiempo, empieza a salir agua con la misma velocidad. En clase se vio que la cantidad
#de sal S(t) que hay en el depósito en el instante t sigue aproximadamente la siguiente
#ecuación diferencial:   
#Se vio también que la solución exacta es S(t) = 60(1 − exp(−t/10)). 
#Use el programa del apartado anterior para resolver el problema con el método de Euler en el intervalo [0, 20]. Repita
#el apartado anterior para este problema.
def ds(t, s):
    """Funcion que define la ecuacion diferencial"""
    return 6- (s/10)
def exactaS(t):
    """Solucion exacta del problema de valor inicial"""
    return 60*(1-exp(-t/10))
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
N = 20 # numero de particiones
s0 = 0. # condicion inicial
figure('Ej1b')
for k in [10, 20, 40, 80, 160]:
    (t, y) = euler(a, b, ds, k, s0) # llamada al metodo de Euler
    ye = exactaS(t) # calculo de la solucion exacta
    # Dibujamos las soluciones
    plot(t, y, '-*') # dibuja la solucion aproximada
    if k == 10:
        leyenda = ['Euler N = 10']
    if k > 10:
        leyenda.append('Euler N = '+ str(k))
    
    if k == 160:
        leyenda.append('Exacta')
        plot(t, ye, 'k') # dibuja la solucion exacta
        xlabel('t')
        ylabel('y')
        legend(leyenda)
        grid(True)
    # Calculo del error cometido
    error = max(abs(y-ye))
    # Resultados
    print('Error',k,': ' + str(error))

    if k!= 160:
        (ts, ys) = euler(a, b, ds, 2*k, s0) # llamada al metodo de Euler
        yes = exactaS(ts) # calculo de la solucion exacta
        errors = max(abs(ys-yes))
        print('Error e',k,'/e',2*k,': ' + str(error/errors))
    print('-----')

show() # muestra la grafica

print('EJERCICIO 1c')
# Suponga que, en el problema anterior, a partir del instante t = 0 en vez de salir agua a razón
#de 2 l/seg lo hace con una velocidad de 3 l/seg, permaneciendo iguales todos los demás
#datos del problema. Modifique la ecuación para tener en cuenta esta variación y aplique el
#método de Euler con N = 2000 a la ecuación resultante en el máximo intervalo de tiempo
#en el que tenga sentido la ecuación. Cuál es, aproximadamente, el máximo de la cantidad
#de sal en el depósito? En qué instante de tiempo se alcanza?
def ds2(t, s):
    """Funcion que define la ecuacion diferencial"""
    return 9- (s/10)
def exactaS2(t):
    """Solucion exacta del problema de valor inicial"""
    return 90*(1-exp(-t/10))
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
N = 20 # numero de particiones
s0 = 0. # condicion inicial
figure('Ej1c')
(t, y) = euler(a, b, ds2, 2000, s0) # llamada al metodo de Euler
ye = exactaS2(t) # calculo de la solucion exacta
# Dibujamos las soluciones
plot(t, y, '-*') # dibuja la solucion aproximada
leyenda.append('Exacta')
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Exacta', 'Euler N = 2000'])
grid(True)
# Calculo del error cometido
error = max(abs(y-ye))
# Resultados
print('Error 2000: ' + str(error))
print('-----')
show() # muestra la grafica

print('EJERCICIO 2a Taylor 2')
# Tomando como modelo el programa anterior, escriba la implementación en Python de los métodos
# de Taylor de orden 2 y 3 para el problema de Cauchy (P). Compruebe el orden repitiendo el
# apartado (a) del ejercicio anterior.
def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 0.5*(t**2 - y)
def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return t**2 - 4*t + 8 - 7.*exp(-0.5*t)

def taylor2(a, b, N, y0):
    """Implementacion del metodo de Taylor2 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Taylor de 2º orden
    for k in range(N):
        dy = 0.5*(t[k]**2 -y[k])
        d2y = t[k]-0.5*dy 
        y[k+1] = y[k]+h*dy + 0.5*h*h*d2y
        t[k+1] = t[k]+h
    
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = 1. # condicion inicial

figure('Ej2a')
for N in [10,20,40,80,160]:
    tini = perf_counter()
    (t, y) = taylor2(a, b, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t) # calculo de la solucion exacta
    plot(t, y, '-*') # dibuja la solucion aproximada
    error = max(abs(y-ye))
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    if N>10:
        cociente=errorold/error
        print('Cociente de errores: ' + str(cociente))
    errorold=error
    print('Paso de malla: ' + str((b-a)/N))
    
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Taylor2, N=10','Taylor2, N=20','Taylor2, N=40','Taylor2, N=80','Taylor2, N=160', 'exacta'])
grid(True)

show() # muestra la grafica
print('EJERCICIO 2b Taylor 3')
def taylor3(a, b, N, y0):
    """Implementacion del metodo de Taylor3 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Taylor de 3er orden
    for k in range(N):
        dy = 0.5*(t[k]**2 -y[k])
        d2y = t[k]-0.5*dy 
        d3y= 1-0.5*d2y
        y[k+1] = y[k]+h*dy + 0.5*h*h*d2y +h*h*h/6 *d3y
        t[k+1] = t[k]+h
    
    return (t, y)


# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = 1. # condicion inicial

figure('Ej2b')
for N in [10,20,40,80,160]:
    tini = perf_counter()
    (t, y) = taylor3(a, b, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t) # calculo de la solucion exacta
    plot(t, y, '-*') # dibuja la solucion aproximada
    error = max(abs(y-ye))
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    if N>10:
        cociente=errorold/error
        print('Cociente de errores: ' + str(cociente))
    errorold=error
    print('Paso de malla: ' + str((b-a)/N))
    
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Taylor3, N=10','Taylor3, N=20','Taylor3, N=40','Taylor3, N=80','Taylor3, N=160', 'exacta'])
grid(True)

show() # muestra la grafica
# Tomando como modelo el programa anterior, escriba la implementación en Python de los métodos
# de Heun, punto medio y RK4. Repita los tres apartados del Ejercicio 1 usando estos métodos
print('EJERCICIO 3 Heun')
def Heun(a, b, fun, N, y0):
    """Implementacion del metodo de Euler en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        t[k+1] = t[k]+h
        ystar = y[k]+h*fun(t[k], y[k])
        y[k+1] = y[k]+0.5*h*(fun(t[k],y[k])+fun(t[k+1],ystar))
        
    
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = 1. # condicion inicial

figure('Ej3a')
for N in [10,20,40,80,160]:
    tini = perf_counter()
    (t, y) = Heun(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t) # calculo de la solucion exacta
    plot(t, y, '-*') # dibuja la solucion aproximada
    error = max(abs(y-ye))
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    if N>10:
        cociente=errorold/error
        print('Cociente de errores: ' + str(cociente))
    errorold=error
    print('Paso de malla: ' + str((b-a)/N))
    
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Heun, N=10','Heun, N=20','Heun, N=40','Heun, N=80','Heun, N=160', 'exacta'])
grid(True)

show() # muestra la grafica

print('EJERCICIO 3 Punto Medio')
def PuntoMedio(a, b, fun, N, y0):
    """Implementacion del metodo de punto medio en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de punto medio
    for k in range(N):
        tstar = t[k]+h*0.5
        ystar = y[k]+0.5*h*fun(t[k], y[k])
        y[k+1] = y[k]+h*fun(tstar,ystar)
        t[k+1] = t[k]+h
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = 1. # condicion inicial

figure('Ej3b')
for N in [10,20,40,80,160]:
    tini = perf_counter()
    (t, y) = PuntoMedio(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t) # calculo de la solucion exacta
    plot(t, y, '-*') # dibuja la solucion aproximada
    error = max(abs(y-ye))
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    if N>10:
        cociente=errorold/error
        print('Cociente de errores: ' + str(cociente))
    errorold=error
    print('Paso de malla: ' + str((b-a)/N))
    
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['PM, N=10','PM, N=20','PM, N=40','PM, N=80','PM, N=160', 'exacta'])
grid(True)

show() # muestra la grafica

print('EJERCICIO 3 RK4')
def RK4(a, b, fun, N, y0):
    """Implementacion del metodo de punto medio en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de punto medio
    for k in range(N):
        k1=fun(t[k],y[k])
        tstar= t[k]+h/2
        k2=fun(tstar,y[k]+0.5*h*k1)
        k3=fun(tstar,y[k]+0.5*h*k2)
        t[k+1]=t[k] +h
        k4=fun(t[k+1],y[k]+h*k3)
        y[k+1]=y[k]+(h/6) *(k1+2*k2+2*k3+k4)
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = 1. # condicion inicial

figure('Ej3c')
for N in [10,20,40,80,160]:
    tini = perf_counter()
    (t, y) = RK4(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t) # calculo de la solucion exacta
    plot(t, y, '-*') # dibuja la solucion aproximada
    error = max(abs(y-ye))
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    if N>10:
        cociente=errorold/error
        print('Cociente de errores: ' + str(cociente))
    errorold=error
    print('Paso de malla: ' + str((b-a)/N))
    
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK4, N=10','RK4, N=20','RK4, N=40','RK4, N=80','RK4, N=160', 'exacta'])
grid(True)

show() # muestra la grafica

print('EJERCICIO 4 EulerSis')
# Escriba implementaciones en Python de los métodos de Euler, Heun, punto medio y RK4 para
# sistemas. Como problema a resolver se considera el modelo depredador/presa de Lotka-Volterra:
# dx/dt = 0.25x − 0.01xy,
# dy/dt = −y + 0.01xy,
# x(0) = 80, y(0) = 30
# Resuelva el problema en el intervalo [0, 20] con N = 20, 40, 80, 160, 320, 640 usando los diferentes
# métodos. Compare en una misma gráfica las trayectorias obtenidas con un mismo método usando
# distintos pasos de tiempo.
# Indicación: Escriba el sistema en forma vectorial y usar la clase array para trabajar con vectores
# y matrices.

def f(t,y):
    f1= 0.25*y[0] - 0.01*y[0]*y[1]
    f2=-y[1] + 0.01*y[0]*y[1]
    return array([f1,f2])
# Datos del problema
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
y0 = array([80,30]) # condicion inicial


def eulerSis(a,b,fun,N,y0): # metodo de euler para sistemas
    """Implementacion del metodo de Euler para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0),N+1])
    t[0] = a
    y[:,0] = y0
    
    for k in range (N):
        y[:,k+1] = y[:,k] + h*fun(t[k],y[:,k])
        t[k+1] = t[k] + h
    return (t,y)

figure('Ej4 EulerSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = eulerSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    
legend(['EulerSis, N=20','EulerSis, N=40','EulerSis, N=80','EulerSis, N=160','EulerSis, N=320','EulerSis, N=640'])
grid(True)
show()

print('EJERCICIO 4 HeunSis')
def heunSis(a, b, fun, N, y0):
    """Implementacion del metodo de Heun para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    ystar = zeros([len(y0),1])
    # Metodo de Heun
    for k in range(N):
        ystar[:,0] = y[:,k] + h*fun(t[k],y[:,k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:,k] + h*1/2*(fun(t[k], y[:,k]) + fun(t[k+1], ystar[:,0]))   
    
    return (t, y)


figure('Ej4 HeunSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = heunSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    
legend(['HeunSis, N=20','HeunSis, N=40','HeunSis, N=80','HeunSis, N=160','HeunSis, N=320','HeunSis, N=640'])
grid(True)
show()

print('EJERCICIO 4 PuntoMedioSis')
def ptoMedioSis(a, b, fun, N, y0):
    """Implementacion del metodo de pto medio para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    ystar = zeros([len(y0),1])
    # Metodo de pto medio
    for k in range(N):
        ystar[:,0] = y[:,k] + h/2*fun(t[k],y[:,k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:,k] + h*fun(t[k]+h/2,ystar[:,0])
        
    return (t, y)

figure('Ej4 PuntoMedioSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = ptoMedioSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    
legend(['PMSis, N=20','PMSis, N=40','PMSis, N=80','PMSis, N=160','PMSis, N=320','PMSis, N=640'])
grid(True)
show()


print('EJERCICIO 4 RK4Sis')
def rk4Sis(a, b, fun, N, y0):
    """Implementacion del metodo de rk4 para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    k1 = zeros([len(y0),1])
    k2 = zeros([len(y0),1])
    k3 = zeros([len(y0),1])
    k4 = zeros([len(y0),1])
    # Metodo de rk4
    for k in range(N):
        k1= fun(t[k],y[:,k])
        k2= fun(t[k]+h/2,y[:,k]+h/2*k1)
        k3= fun(t[k]+h/2,y[:,k]+h/2*k2)
        k4= fun(t[k]+h,y[:,k]+h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:,k] + h/6*(k1 +2*k2 + 2*k3 + k4)
        
    return (t, y)

figure('Ej4 RK4Sis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    
legend(['RK4Sis, N=20','RK4Sis, N=40','RK4Sis, N=80','RK4Sis, N=160','RK4Sis, N=320','RK4Sis, N=640'])
grid(True)
show()

print('EJERCICIO 5 EulerSis')
# Se considera el problema de Cauchy:
# x''+ 20x' + 101x = 0,
# x(0) = 1, x'(0)= −10.
# cuya solución es x(t) = exp(-10t)cos(t). 
# Reescriba la ecuación como un sistema de dos ecuaciones de primer orden y aplique los métodos del 
# ejercicio anterior en el intervalo [0, 7] con N = 20, 40, 80, 160, 320, 640. Compare en una misma gráfica
# la solución exacta x(t) con las aproximaciones obtenidas con un mismo método usando distintos pasos de tiempo. 
# Se observan resultados anómalos?
def f(t, y):
    f1 = y[1]
    f2 = -20*y[1] - 101*y[0]
    return array([f1, f2])
def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return exp(-10*t)*cos(t)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 7. # extremo superior del intervalo
y0 = array([1,-10]) # condicion inicial

figure('Ej5 EulerSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = eulerSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    

legend(['EulerSis, N=20','EulerSis, N=40','EulerSis, N=80','EulerSis, N=160','EulerSis, N=320','EulerSis, N=640'])
plot(t, exacta(t), 'k') # dibuja la solucion exacta
grid(True)
show()

print('EJERCICIO 5 HeunSis')

figure('Ej5 HeunSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = heunSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    

legend(['HeunSis, N=20','HeunSis, N=40','HeunSis, N=80','HeunSis, N=160','HeunSis, N=320','HeunSis, N=640'])
plot(t, exacta(t), 'k') # dibuja la solucion exacta
grid(True)
show()

print('EJERCICIO 5 PuntoMedioSis')

figure('Ej5 PuntoMedioSis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = ptoMedioSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    
legend(['PMSis, N=20','PMSis, N=40','PMSis, N=80','PMSis, N=160','PMSis, N=320','PMSis, N=640'])
plot(t, exacta(t), 'k') # dibuja la solucion exacta
grid(True)
show()


print('EJERCICIO 5 RK4Sis')

figure('Ej5 RK4Sis')
for N in [20, 40, 80, 160, 320, 640]:
    tini = perf_counter()
    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    subplot(211)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    subplot(212)
    plot(y[0,:],y[1,:]) #grafica soluciones x,y
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    

legend(['RK4Sis, N=20','RK4Sis, N=40','RK4Sis, N=80','RK4Sis, N=160','RK4Sis, N=320','RK4Sis, N=640'])
plot(t, exacta(t), 'k') # dibuja la solucion exacta
grid(True)
show()

print('EJERCICIO 6')
# Modelo SIR Para simular la evolución de una epidemia, se consideran las siguientes tres clases disjuntas de individuos:
# S: individuos susceptibles de ser contagiados.
# I: individuos infectados.
# R: individuos que no están en ninguna de las dos clases anteriores, ya sea por aislamiento, recuperación (y, consecuentemente, inmunidad), o por muerte.
# Sean x(t), y(t) y z(t) el número de individuos (medido en millares) que están en el instante t (medido en días) en las clases S, I y R,
 # respectivamente. Se considera el siguiente modelo para simular la evolución de la epidemia:
# dx/dt = -k1xy - k3x
# dy/dt = k1xy - k2y
# dz/dt = k3x + k2y
# siendo k1 la tasa de contagios, k2 la tasa de eliminación de infectados por las diversas causas posibles y k3 la tasa de mortalidad de los individuos sanos. 
# a) Supongamos que el número total de individuos en el instante t = 0 es N = 100 y que en el instante inicial hay 1 millar 
# contagiados y 99 millares de susceptibles. Supongamos además que k1 = 2 × 10−3, k2 = 5 · 10−3 y k3 = 10−5. 
# Aproxime la solución del sistema de ecuaciones diferenciales en el intervalo de tiempo [0, 60] usando el método RK4 con paso h = 0.05.
# Represente en una sola gráfica la evolución de x(t), y(t), z(t). Comente brevemente los resultados que parecen deducirse sobre la evolución de la epidemia
# b) Repita el apartado anterior con todos los datos iguales, pero con k2 = 5 · 10−2: compare las diferencias.
# c) Modifique el modelo para que contemple la posibilidad de que individuos infectados puedan pasar a ser de nuevo susceptibles por pérdida de inmunidad, siendo la velocidad de
# este proceso proporcional al número de infectados: k4 y. Simule la evolución con los datos del apartado anterior tomando k4 = 0.5 k2: compare los resultados.

c1 = 2.e-3   ### c1 tasa de contagios
c2 = 5.e-3   ### c2 tasa de eliminacion de infectados
c3 = 1.e-3     ### c3 tasa de mortalidad
def f(t,y):
    f1 = -c1*y[0]*y[1] - c3*y[0]
    f2 = c1*y[0]*y[1] - c2*y[1]
    f3 = c3*y[0] + c2*y[1]
    return array([f1,f2,f3])

a = 0
b = 60
y0 = array([99,1,0])
h=0.05
N= 1200

figure('Ej5a')
tini = perf_counter()
(t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de rk4 para sistemas
tfin=perf_counter()
plot(t, y[0], t, y[1], t, y[2]) 
legend(['Individuos suceptibles','Individuos infectados','Individuos recuperados'])

print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((b-a)/N))
show()

### Apartado b
c2 = 5.e-2   ### c2 tasa de eliminacion de infectados

figure('Ej5b')
tini = perf_counter()
(t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de rk4 para sistemas
tfin=perf_counter()
plot(t, y[0], t, y[1], t, y[2]) 
legend(['Individuos suceptibles','Individuos infectados','Individuos recuperados'])

print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((b-a)/N))
show()

### Apartado c
c2= 5.e-3
c4= 5.e-2

def f(t,y):
    f1 = -c1*y[0]*y[1] - c3*y[0] + c4*y[1]
    f2 = c1*y[0]*y[1] - c2*y[1] -c4*y[1]
    f3 = c3*y[0] + c2*y[1]
    return array([f1,f2,f3])

figure('Ej5c')
tini = perf_counter()
(t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de rk4 para sistemas
tfin=perf_counter()
plot(t, y[0], t, y[1], t, y[2]) 
legend(['Individuos suceptibles','Individuos infectados','Individuos recuperados'])

print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((b-a)/N))
show()

print('EJERCICIO 7a')
# Modelo de cohete propulsado. El sigiuiente problema de Cauchy modela el movimiento de un cohete que se lanza verticalmente:
# dz/dt = v
# dv/dt = -g + T/(M + mf) - C*v*|v|/(M + mf) + α*T*v/(M + mf) 
# dm/dt = -αT
# z(0) = 0
# v(0) = v0
# mf(0) = mf,0
# siendo z(t) la altura del centro de gravedad del cohete en el instante t (se toma como altura 0 la posición inicial); 
# v(t) la velocidad del cohete; T(t) la fuerza ejercida por el motor; mf(t) la masa de combustible; M la masa del cohete con el depósito
# de combustible vacío; C el coeficiente de fricción con el aire; v0 la velocidad de lanzamient y mf ,0 la masa inicial de combustible.
# La velocidad a la que se reduce la masa de combustible se supone proporcional a la fuerza que ejerce el motor, siendo α la constante
# de proporcionalidad. Supondremos los siguientes valores: g = 9.81, M = 7.5, α = 0.02. Supondremos además que la fuerza que ejerce el
# cohete es una función de la forma: T(t) = T0 si mf > 0,0 si mf = 0, siendo T0 un valor constante.
# a) Esciriba un programa Python que permita resolver el problema planteado usando el método RK4 con paso h mientras se cumpla z > 0.
# b) Aplique el programa para resolver el problema con condición inicial v0 = 50, mf ,0 = 7.5, paso h = 0.05 y las siguientes elecciones
#     de T0 y C: (i) T0 = C = 0; (ii) T = 0, C = 0.02; (iii) T0 = 50, C = 0.02
# c) Calcule la altura máxima alcanzada por el cohete en cada caso y estime el tiempo de caída. Compare las gráficas de la altura del
#     cohete frente al tiempo obtenida en los tres casos. En el caso (c) dibuje también la gráfica de la masa de combustible frente
#     al tiempo y calcule en qué momento se acaba el combustible.

def cohete_propulsado(T0, C):
    # Condiciones iniciales
    z0 = 0
    mf0 = 7.5
    v0 = 50  
    y0 = array([z0, v0, mf0])
    
    def cohete(t, y):
        g = 9.81
        M = 7.5
        alpha = 0.02
        # T0 = T0 
        # C = C 
        z, v, mf = y
        
        if mf > 0:
            T = T0
        else:
            T = 1e-9 # evitar dividir por 0
    
        # limitar la velocidad para evitar sobrecarga
        v_max = 3000
        v = max(min(v, v_max), -v_max)
            
        dzdt = v
        dvdt = -g + T/(M + mf) - C*v*abs(v)/(M + mf) + alpha*T*v/(M + mf)
        dmdt = -alpha*T
        
        return array([dzdt, dvdt, dmdt])
    
    # Intervalo y número de particiones
    a = 0
    b = 50
    N = 1000
    h = (b-a)/N # Paso de malla = 0.05
    print('Paso de malla usado:' + str(h))
    # Resolución mediante RK4
    t, y = rk4Sis(a, b, cohete, N, y0)

    # Extracción de las soluciones
    z = y[0, :]
    v = y[1, :]
    mf = y[2, :]

    # Representación gráfica
    plt.plot(t, z)
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Altura (m)')
    plt.show()
    
    return t, y

print('EJERCICIO 7b')
figure('Ej7b')
# Caso (i) T0 = C = 0
t1, y1 = cohete_propulsado(T0=0, C=0)

# Caso (ii) T = 0, C = 0.02
t2, y2 = cohete_propulsado(T0=0, C=0.02)

# Caso (iii) T0 = 50, C = 0.02
t3, y3 = cohete_propulsado(T0=50, C=0.02)
fig, ax = plt.subplots()

ax.plot(t1, y1[0], label='T0=0, C=0')
ax.plot(t2, y2[0], label='T0=0, C=0.02')
ax.plot(t3, y3[0], label='T0=50, C=0.02')

ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Altura (m)')
ax.legend()

show()

print('EJERCICIO 7c')
def cohete_propulsado(T0, C):
    # Condiciones iniciales
    z0 = 0
    mf0 = 7.5
    v0 = 50  
    y0 = array([z0, v0, mf0])
    
    def cohete(t, y):
        g = 9.81
        M = 7.5
        alpha = 0.02
        # T0 = T0 
        # C = C 
        z, v, mf = y
        
        if mf > 0:
            T = T0
        else:
            T = 1e-9 # evitar dividir por 0
    
        # limitar la velocidad para evitar sobrecarga
        v_max = 3000
        v = max(min(v, v_max), -v_max)
            
        dzdt = v
        dvdt = -g + T/(M + mf) - C*v*abs(v)/(M + mf) + alpha*T*v/(M + mf)
        dmdt = -alpha*T
        
        return array([dzdt, dvdt, dmdt])
    
    # Intervalo y número de particiones
    a = 0
    b = 50
    N = 1000
    h = (b-a)/N # Paso de malla = 0.05
    print('Paso de malla usado:' + str(h))
    # Resolución mediante RK4
    t, y = rk4Sis(a, b, cohete, N, y0)

    # Extracción de las soluciones
    z = y[0, :]
    v = y[1, :]
    mf = y[2, :]

    # Representación gráfica de la altura
    plt.subplot(2, 1, 1)
    plt.plot(t, z)
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Altura (m)')

    # Representación gráfica de la velocidad
    plt.subplot(2, 1, 2)
    plt.plot(t, v)
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Velocidad (m/s)')
    
    plt.show()
    
    return t, y, v

print('EJERCICIO 7c')
figure('Ej7c')
# Caso (i) T0 = C = 0
t1, y1, v1 = cohete_propulsado(T0=0, C=0)

# Caso (ii) T = 0, C = 0.02
t2, y2, v2 = cohete_propulsado(T0=0, C=0.02)

# Caso (iii) T0 = 50, C = 0.02
t3, y3, v3 = cohete_propulsado(T0=50, C=0.02)

fig, axs = plt.subplots(2)

axs[0].plot(t1, y1[0], label='T0=0, C=0')
axs[0].plot(t2, y2[0], label='T0=0, C=0.02')








