# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pylab import *
from time import perf_counter

# Mide el tiempo CPU que tarda en metodos numericos



print('EJERCICIO 1')
# Para comprobar que las implementaciones de los métodos de Heun, punto medio y RK4 son
# correctas, ejecutad los códigos para resolver el problema
# y' = −y + 2 sen(t), t ∈ [0, 10],
# y(0) = π,
# cuya solución exacta es la función y(t) = (π + 1)e^(−t) + sen(t) − cos(t), con N = 50 particiones.
# Los errores que deberían obtenerse son los siguientes:
# Euler: 0.193460359895
# Heun: 0.0120152817472
# Punto medio: 0.0146795451937
# RK4: 2.4240784029 × 10−5
# Indicación: El error se calcula como max||y(tk) − yk||   0<=t<=N
def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return -y + 2*sin(t)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return (pi + 1)*exp(-t) + sin(t) - cos(t)


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
N = 50 # numero de particiones
y0 = pi # condicion inicial

tini = perf_counter() # Tiempo CPU antes de llamar al metodo

(t, y) = euler(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter() # Tiempo CPU despues de llamar al metodo

ye = exacta(t) # calculo de la solucion exacta

print('Euler')
# Dibujamos las soluciones
figure('Ej1E')
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
tini = perf_counter() # Tiempo CPU antes de llamar al metodo

(t, y) = Heun(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter() # Tiempo CPU despues de llamar al metodo

ye = exacta(t) # calculo de la solucion exacta
print('Heun')
# Dibujamos las soluciones
figure('Ej1H')
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Heun', 'exacta'])
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

tini = perf_counter() # Tiempo CPU antes de llamar al metodo

(t, y) = PuntoMedio(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter() # Tiempo CPU despues de llamar al metodo

ye = exacta(t) # calculo de la solucion exacta
print('PuntoMedio')
# Dibujamos las soluciones
figure('Ej1PM')
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['PM', 'exacta'])
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

tini = perf_counter() # Tiempo CPU antes de llamar al metodo

(t, y) = RK4(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter() # Tiempo CPU despues de llamar al metodo

ye = exacta(t) # calculo de la solucion exacta

print('PuntoMedio')
# Dibujamos las soluciones
figure('Ej1RK')
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK', 'exacta'])
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

print('EJERCICIO 2')
# Para hacer lo propio con las implementaciones para sistemas, considerad el problema
# x'(t) = 3x(t) − 2y(t),
# y'(t) = −x(t) + 3y(t) − 2z(t),
# z'(t) = −y(t) + 3z(t),
# en el intervalo [0, 1] con condiciones iniciales
# x(0) = 1, y(0) = 0, z(0) = −1,
# y resolvedlo tomando N = 50 particiones. La solución exacta es
# x(t) = −1/4 *e^(5t) + 3/2 *e(3t) − 1/4 *e^t,
# y(t) = 1/4 *e^(5t) − 1/4 *e^(t),
# z(t) = −1/8 e^(5t) − 3/4 *e^(3t) − 1/8 *e^(t).
# Los errores en (x, y, z) que deberían obtenerse son los siguientes:
# Euler:
# • Error x: 5.26417603327,
# • Error y: 7.74890310587,
# • Error: z: 5.13016198893.
# Heun:
# • Error x: 0.234062642044,
# • Error y: 0.28577784976,
# • Error: z: 0.168835786932.
# Punto medio:
# • Error x: 0.234062642044,
# • Error y: 0.28577784976,
# • Error: z: 0.168835786932.
# RK4:
# • Error x: 0.000132965132994,
# • Error y: 0.000142249175603,
# • Error: z: 7.57683913406 × 10−5

# Indicación: El error se calcula como en el apartado anterior, componente a componente

def f(t,y):
    f1= 3*y[0] - 2*y[1]
    f2= - y[0] + 3*y[1] - 2*y[2]
    f3= - y[1] + 3*y[2]
    return array([f1,f2, f3])

def exacta(t):
    f1= -exp(5*t)/4. + 3*exp(3*t)/2. - exp(t)/4.
    f2= exp(5*t)/4. - exp(t)/4.
    f3= -exp(5*t)/8. - 3*exp(3*t)/4. - exp(t)/8.
    return(array([f1,f2, f3]))


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

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 1. # extremo superior del intervalo
N = 50 # numero de particiones
y0 = [1, 0, -1] # condicion inicial

figure('Ej2EuLSis')

tini = perf_counter()
(t, y) = eulerSis(a, b, f, N, y0) # llamada al metodo de Euler
tfin=perf_counter()

ye = exacta(t) # calculo de la solucion exacta

errorx = max(abs(y[0]-ye[0]))
errory = max(abs(y[1]-ye[1]))
errorz = max(abs(y[2]-ye[2]))

subplot(211)
grid(True)
plot(t, y[0], t , y[1], t, y[2]) #grafica t
subplot(212)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(y[0,:], y[1,:], y[2,:])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# plot(y[0,:],y[1,:], y[2, :]) #grafica soluciones x,y
# legend(['x', 'y', 'z'])
# Graficar las soluciones
# plt.plot(t, y[0,:], label='x')
# plt.plot(t, y[1,:], label='y')
# plt.plot(t, y[2,:], label='z')
# plt.legend()
# plt.show()


print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error x: ' + str(errorx))
print('Error y: ' + str(errory))
print('Error z: ' + str(errorz))
print('Paso de malla: ' + str((b-a)/N))


grid(True)
show()


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

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 1. # extremo superior del intervalo
N = 50 # numero de particiones
y0 = [1, 0, -1] # condicion inicial

figure('Ej2HeunSis')

tini = perf_counter()
(t, y) = heunSis(a, b, f, N, y0) # llamada al metodo de Euler
tfin=perf_counter()

ye = exacta(t) # calculo de la solucion exacta

errorx = max(abs(y[0]-ye[0]))
errory = max(abs(y[1]-ye[1]))
errorz = max(abs(y[2]-ye[2]))

subplot(211)
grid(True)
plot(t, y[0], t , y[1], t, y[2]) #grafica t
subplot(212)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(y[0,:], y[1,:], y[2,:])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# plot(y[0,:],y[1,:], y[2, :]) #grafica soluciones x,y
# legend(['x', 'y', 'z'])
# Graficar las soluciones
# plt.plot(t, y[0,:], label='x')
# plt.plot(t, y[1,:], label='y')
# plt.plot(t, y[2,:], label='z')
# plt.legend()
# plt.show()


print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error x: ' + str(errorx))
print('Error y: ' + str(errory))
print('Error z: ' + str(errorz))
print('Paso de malla: ' + str((b-a)/N))


grid(True)
show()



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

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 1. # extremo superior del intervalo
N = 50 # numero de particiones
y0 = [1, 0, -1] # condicion inicial

figure('Ej2PMSis')

tini = perf_counter()
(t, y) = ptoMedioSis(a, b, f, N, y0) # llamada al metodo de Euler
tfin=perf_counter()

ye = exacta(t) # calculo de la solucion exacta

errorx = max(abs(y[0]-ye[0]))
errory = max(abs(y[1]-ye[1]))
errorz = max(abs(y[2]-ye[2]))

subplot(211)
grid(True)
plot(t, y[0], t , y[1], t, y[2]) #grafica t
subplot(212)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(y[0,:], y[1,:], y[2,:])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# plot(y[0,:],y[1,:], y[2, :]) #grafica soluciones x,y
# legend(['x', 'y', 'z'])
# Graficar las soluciones
# plt.plot(t, y[0,:], label='x')
# plt.plot(t, y[1,:], label='y')
# plt.plot(t, y[2,:], label='z')
# plt.legend()
# plt.show()


print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error x: ' + str(errorx))
print('Error y: ' + str(errory))
print('Error z: ' + str(errorz))
print('Paso de malla: ' + str((b-a)/N))


grid(True)
show()

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


# Datos del problema
a = 0. # extremo inferior del intervalo
b = 1. # extremo superior del intervalo
N = 50 # numero de particiones
y0 = [1, 0, -1] # condicion inicial

figure('Ej2RK4Sis')

tini = perf_counter()
(t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de Euler
tfin=perf_counter()

ye = exacta(t) # calculo de la solucion exacta

errorx = max(abs(y[0]-ye[0]))
errory = max(abs(y[1]-ye[1]))
errorz = max(abs(y[2]-ye[2]))

subplot(211)
grid(True)
plot(t, y[0], t , y[1], t, y[2]) #grafica t
subplot(212)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(y[0,:], y[1,:], y[2,:])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# plot(y[0,:],y[1,:], y[2, :]) #grafica soluciones x,y
# legend(['x', 'y', 'z'])
# Graficar las soluciones
# plt.plot(t, y[0,:], label='x')
# plt.plot(t, y[1,:], label='y')
# plt.plot(t, y[2,:], label='z')
# plt.legend()
# plt.show()


print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error x: ' + str(errorx))
print('Error y: ' + str(errory))
print('Error z: ' + str(errorz))
print('Paso de malla: ' + str((b-a)/N))


grid(True)
show()































