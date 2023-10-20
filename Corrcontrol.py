from pylab import *
from time import perf_counter

# Isidro Javier Garcia Fernandez
# Doble grado matematicas e ingenieria informatica

# Un cohete de M = 12kg de peso cae a tierra siguiendo una trayectoria vertical.
# En el instante t = 0s se encuentra a una altura h = 1000m y cae a una velocidad de 25m/s.
# En ese momento se enciende un cohete retropropulsor que ejerce una fuerza hacia arriba de intensidad T = 119 (se desprecian las variaciones
# de masa debidas al consumo de combustible). Se desea aproximar la velocidad del cohete al llegar al suelo. Para ello, se considera el problema de Cauchy:
# z''(t) = -g + T/m - C/m*z'(t)*|z'(t)|
# z(0) = 1000, z'(0) = -25
# donde z(t) representa la altura del cohete en el instante t, g = 9.81 /s^2 es la aceleracion de la gravedad y C = 0.01 es el coeficiente de friccion.

# 1. Escriba un programa en python que resuelva el problema de Cauchy usando el método RK4 con paso h= 0.1. Se resolverán las ecuaciones
# mientras que la altura del cohete sea mayor que 0, o por seguridad, hasta que se alcance un tiempo máximo de 1000s.
# Si el programa se detiene porque se ha alcanzado el tiempo máximo, tendrá que advertirse en pantalla.
# Dibuje, en una única ventana gráfica (usando el comando subplot) las gráficas de la altura y de la velocidad del cohete frente al tiempo.
# Ponga en pantalla el tiempo de calculo, el tiempo que ha tardado el cohete en llegar a tierra y la velocidad del cohete en el último instante de tiempo calculado.

g = 9.81
M = 12
T = 119
C = 0.01
maxstep = 10000 #el programa se detiene si se han dado 10000 pasos para evitar bucles infinitos

def fun(t,z):
    #z[0]:= altura; z[1]:=velocidad; 
    
    f0 = z[1]
    f1 = -g + T/M - C/M*z[1]*abs(z[1])
    
    return array([f0,f1])

def RK4Sis(a, fun, y0, h):
    """Implementacion del metodo de rk4 para el problema del cohete
    usando paso h y condicion inicial y0"""
    k = 0  # contador de pasos
    t = array([a]) # inicializacion del vector de nodos CORRECCION
    y = zeros([len(y0), 1])  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[:, 0] = y0  # valor inicial
    k1 = zeros([len(y0), 1])
    k2 = zeros([len(y0), 1])
    k3 = zeros([len(y0), 1])
    k4 = zeros([len(y0), 1])
    

    # Condiciones para el bucle
    tini = perf_counter()   # CORRECCIÓN
    tfin = tini             # CORRECCIÓN

    # Metodo de rk4
    while (tfin - tini <= 1000) and y[0,k] > 0:
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        
        t = append(t, t[k]+h) # CORRECCIÓN
        
        t[k+1] = t[k]+h
        y = column_stack((y, y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4))) # CORRECCION
        k = k+1
        
        tfin = perf_counter() # CORRECCIÓN

    if (tfin - tini >= 1000):
        print("El programa se ha parado porque hemos alcanzado 1000s")
        
    k -= 1
    return (t, y, k) # CORRECCION (devolver k)

# Condiciones iniciales
h = 0.1 # paso de malla
a = 0 # tiempo inicial
z0 = array([1000, -25]) # altura inicial y velocidad inicial


tinicial = perf_counter()
t, z, k = RK4Sis(a, fun, z0, h) 
tfinal = perf_counter()
print('EJERCICIO 1')
print('------------')
print("Tiempo de calculo: ", abs(tfinal-tinicial))
print("Tiempo que ha tardado el cohete en llegar a tierra: ", t[k]) # CORRECCION
print("Velocidad del cohete en el ultimo instante de tiempo calculado: ", z[1, k]) # CORRECCION
print('------------')
# Dibujamos las graficas
figure(1)
subplot(2,1,1)
plot(t,z[0])
xlabel("Tiempo")
ylabel("Altura")
subplot(2,1,2)
plot(t,z[1])
xlabel("Tiempo")
ylabel("Velocidad")
show()



# Ejercicio 2
def rk45(a, fun, y0, h0, tol): #de quinto orden necesita al menos 6 etapas
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""
    hmin = 1e-5 # paso de malla minimo
    hmax = 0.05 # paso de malla maximo

    
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
    y = y0.reshape(len(y0), 1) # soluciones
    h = array([h0]) # pasos de malla
    K = zeros((len(y0), q)) #array que contiene k1, k2, k3 que usa Runge-Kutta para avanzar
    k = 0 # contador de iteraciones
    
    # Condiciones para el bucle
    tini = perf_counter()       #CORRECCION
    tfin = tini                 #CORRECCION
    
    
    while (tfin - tini <= 1000) and y[0,k] > 0:  # CORRECCION
    
        for i in range(q): #calcular K1, K2, K3
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]* dot(A[i, :], transpose(K)))
        incrlow = dot(B, transpose(K)) # metodo de orden 2 #phi
        # B es b1, b2, b3, K es k1, k2, k3 multiplica elemento a elemento y los suma
        incrhigh = dot(BB, transpose(K)) # metodo de orden 3 #phiEstrella
        # BB es b1*, b2*, b3*, K es k1, k2, k3 multiplica elemento a elemento y los suma
        
#       error = h[k]*(incrhigh-incrlow) # estimacion del error
#       y = append(y, y[k]+h[k]*incrlow) # y_(k+1)
#       t = append(t, t[k]+h[k]); # t_(k+1)
#       hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1)) # h_(k+1)
#       hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax

        
        error = norm(h[k]*(incrhigh-incrlow), inf) # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow)) # y_(k+1)
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1)) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        
        
        
        h = append(h, hnew)
        k += 1
        
        tfin = perf_counter() # CORRECCION
        
        
    if (tfin - tini >= 1000):
        print("El programa se ha parado porque hemos alcanzado 1000s")

    # print("Numero de iteraciones: " + str(k))
    return (t, y, h, k)

# Condiciones iniciales
h0 = 0.01 # paso de malla
a = 0 # tiempo inicial
tol = 1e-8 # tolerancia de error
z0 = array([1000, -25]) # altura inicial y velocidad inicial

tinicial = perf_counter()
(t, z, h, k) = rk45(a, fun, z0, h0, tol) # llamada al metodo RK2(3)
tfinal = perf_counter()
print('EJERCICIO 2')
print('------------')
print("Tiempo de calculo: ", abs(tfinal-tinicial))
print("Tiempo que ha tardado el cohete en llegar a tierra: ", t[k])
print("Velocidad del cohete en el ultimo instante de tiempo calculado: ", z[1, k])

# Dibujamos las graficas
figure(2)
subplot(311)
plot(t,z[0])
xlabel("Tiempo")
ylabel("Altura")
subplot(312)
plot(t,z[1])
xlabel("Tiempo")
ylabel("Velocidad")
subplot(313)
plot(t, h)# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
show()

hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
