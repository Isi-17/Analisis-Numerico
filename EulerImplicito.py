from pylab import *
from time import perf_counter

# x'' + 20x' + 101x = 0  t€ [0, 7]
# x(0) = 1, x'(0) = -10
# x(t) = exp(-10*t)*cos(t)

# 1. Escriba el problema en forma equivalente como un sistema de dos ecuaciones de primer orden.
# En los siguientes apartados, basta dibujar las gráficas y calcular los errores para la variable x(t)

# 2. Aplique el método RK4 con pasos h = 7/N, N = 25, 100, 400. Determine el error cometido.
# Dibuje las gráficas de las soluciones numéricas junto con la de la solución exacta. Extraiga conclusiones.

# 3. Idem para el método de Euler implicito. Observe que, en cada etapa del método, hay quye resolver un 
# sistema lineal de dos ecuaciones y dos incógnitas: calcule a mano la solución exacta de dicho sistema y 
# use su expresión para implementar el método.

#  Ejercicio 1
def f(t, y):
    return array([y[1], -20*y[1] + -101*y[0]])

def exacta(t):
    return exp(-10*t)*cos(t)

#Datos del problema
a = 0
b = 7
y0 = array([1, -10])
malla = [25, 100, 400]


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

figure('Ej1 RK4Sis')
for N in malla:
    tini = perf_counter()
    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t)
    error = abs(y[0,:] - ye)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    print('-----')
    print('N = ' + str(N))
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('Error maximo: ' + str(max(error)))
    
legend(['RK4Sis, N=25','RK4Sis, N=100','RK4Sis, N=400'])
show()

def eulerImplicitoSis(a,b,fun,N,y0):
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    
    for k in range(N):
        y[1,k+1] = (y[1,k] -101*h*y[0,k])/(1+20*h+101*h*h)
        y[0,k+1] = y[0,k] +h*y[1,k+1]
        t[k+1] = t[k] + h
    
    return(t,y)

figure('Ej1 EulerImplicitoSis')
for N in malla:
    tini = perf_counter()
    (t, y) = eulerImplicitoSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye = exacta(t)
    error = abs(y[0,:] - ye)
    grid(True)
    plot(t, y[0], t , y[1]) #grafica t
    print('-----')
    print('N = ' + str(N))
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('Error maximo: ' + str(max(error)))

legend(['EulerImplicitoSis, N=25','EulerImplicitoSis, N=100','EulerImplicitoSis, N=400'])
show()



