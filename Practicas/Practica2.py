# Resolucion del problema de valor inicial
# y'=f(t,y), y(t0)=y0,
# mediante el metodo RK2(3).

from pylab import *
from time import perf_counter

def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return -y+exp(-t)*cos(t)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return exp(-t)*sin(t)

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
    
    # print("Numero de iteraciones: " + str(k))
    return (t, y, h)


print('EJERCICIO 1')
    
# Datos del problema
a = 0. # extremo inferior del intervalo
b = 30. # extremo superior del intervalo
y0 = 0. # condicion inicial
h0 = 0.1 #paso inicial
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
    
    # print("Numero de iteraciones: " + str(k))
    return (t, y, h)

print('EJERCICIO 2')

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

print('EJERCICIO 3')
# Modifique el programa correspondiente el método RK2(3) a fin de poderlo aplicar a la resolución
# de sistemas. (Observación: en el caso de un sistema tanto el error de discretización local εk como
# su estimador eεk son vectores: sustituya |eεk| por |eεk||∞ en la expresión del nuevo paso hk+1.)
# Aplique el programa a la resolución de los problemas 4, 5, 6(a), 7(b) y 7(c) de la Práctica 1 con
# tolerancia 10−4 y 10−6. Además de las gráficas de las soluciones, dibuje en cada caso otra que
# muestre la evolución del paso de malla.
# dx/dt = 0.25x − 0.01xy,
# dy/dt = −y + 0.01xy,
# x(0) = 80, y(0) = 30
# Intervalo [0, 20]
# Ejercicio 4 Práctica 1
def f(t,y):
    f1= 0.25*y[0] - 0.01*y[0]*y[1]
    f2=-y[1] + 0.01*y[0]*y[1]
    return array([f1,f2])

def exacta(t):
    f1= 80 *exp(t/4)
    f2= 250/(1 + sqrt(64001)*exp(t/8))
    return(array([f1,f2]))

# x(t) = 80e^(0.25t)
# y(t) = 1250/[1 + sqrt(64001)e^(0.125t)]

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
    # print("Numero de iteraciones: " + str(k))   
    return (t, y, h)




# Datos del problema
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
y0 = array([80,30]) # condicion inicial
# y0 = y0.reshape(2,1)
h0 = (b-a)/2000 #paso inicial
tol = 1.e-6 #tolerancia

tini = perf_counter()
(t, y, h) = rkSistemas23(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()

figure('Ej3a')
# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 


subplot(311)
plot(t,y[0],t,y[1]) 
subplot(312)
plot(y[0],y[1])
subplot(313)
plot(t, h) 

print('-----')
# Calculo del error cometido
error1= max(abs(y[0]-exacta(t)[0]))
error2= max(abs(y[1]-exacta(t)[1]))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')

print("No. de iteraciones: " + str(len(y[0])-1))
print(error1)
print(error2)
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')

# def f(t,y):
#     f1= 0.25*y[0] - 0.01*y[0]*y[1]
#     f2=-y[1] + 0.01*y[0]*y[1]
#     f3 = 0
#     return array([f1,f2, f3])

# def exacta(t):
#     f1= 80 *exp(t/4)
#     f2= 250/(1 + sqrt(64001)*exp(t/8))
#     return(array([f1,f2]))

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
    # print("Numero de iteraciones: " + str(k))        
    return (t, y, h)

# # Datos del problema
# a = 0. # extremo inferior del intervalo
# b = 1. # extremo superior del intervalo
# y0 = array([1,0,-1]) # condicion inicial
# # y0 = y0.reshape(3,1)
# h0 = (b-a)/2000 #paso inicial
# tol = 1.e-6 #tolerancia

# tini = perf_counter()
# (t, y, h) = rkSistemas45(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
# tfin = perf_counter()

# figure('Ej3b')
# # calculo de la solucion exacta
# te = linspace(a,b,200)
# ye = exacta(te) 
     
# subplot(311)
# plot(t,y[0],t,y[1],t,y[2]) 
# subplot(312)
# plot(y[0],y[1])
# subplot(313)
# plot(t, h)
# print('-----')
# # Calculo del error cometido
# error1= max(abs(y[0]-exacta(t)[0]))
# error2= max(abs(y[1]-exacta(t)[1]))
# error3= max(abs(y[2]-exacta(t)[2]))
# hn = min(h[:-2]) # minimo valor utilizado del paso de malla
# hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# # se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# # Resultados
# print('-----')

# print("No. de iteraciones: " + str(len(y[0])-1))
# print(error1)
# print(error2)
# print(error3)
# print('Tiempo CPU: ' + str(tfin-tini))
# print("Paso de malla minimo: " + str(hn))
# print("Paso de malla maximo: " + str(hm))
# print(sum(h))
# print('-----')






