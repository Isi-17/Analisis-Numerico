from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *
from time import perf_counter
from numpy.linalg import eig, norm
# Se considera el problema de Cauchy:

# x' = y,
# y' = −αx − βy
# x(0) = 1,
# y(0) = 0.

# 1. Resuelva el problema con las siguientes elecciones de parámetros
# a) α = 10, β = 1;
# b) α = 15, β = 3.
# en el intervalo [0, 20]. Para ello, use el método el método Adams-Bashforth de 4 pasos (AB4) 
# yk+1 = yk + h/24 (55 fk − 59 fk−1 + 37 fk−2 − 9 fk−3)   , k = 3, 4, . . .
# con una partición uniforme de N = 200 subintervalos. Use el método RK4 para arrancar. Para cada una de las elecciones de parámetros dibuje:
#  -> la gráfica de x frente a t, la gráfica de y frente a t y la trayectoria en el plano de fases (deben de ir separadas pero en una única figura:
#     use para ello el comando subplot)
#  -> la frontera de la región de estabilidad absoluta del método en el plano complejo y los puntos donde se sitúan hλi , i = 1, 2 en dicho plano,
#     siendo h = 20/N y λi, i = 1, 2 los autovalores de la matriz del sistema.

# Si cree que la posición de hλi , i = 1, 2 en relación a la frontera de la región de estabilidad explica de alguna manera los resultados obtenidos, coméntelo

a = 15
b = 3

# y[0] = t
# y[1] = x
# y[2] = y

def fun2(t, y):
    f1 = y[1]
    f2 = -a*y[0] - b*y[1]
    return array([f1, f2])

# def exacta2(t):
#     f1 = 1/2 * exp(-t) * (exp((sqrt(39) - 1) * t) + exp(-(sqrt(39) + 1) * t))
#     f2 = 1/2 * exp(-t) * (sqrt(39) * exp((sqrt(39) - 1) * t) - sqrt(39) * exp(-(sqrt(39) + 1) * t))
#     return array([f1, f2])

ini = 0
fin = 20
N = 200
h = (fin - ini) / N
y0 = array([1, 0])

def AB4Sis(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])

    ## Metodo de Runge-Kutta de orden 4 Sistemas
    for k in range(3):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k] + h/2, y[:, k] + h/2 *k1)
        k3 = fun(t[k] + h/2*h, y[:, k] + h/2 *k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
    
    for k in range(3, N):
        y[:, k+1] = y[:, k] + h/24 * (55 * f[:, k] - 59 * f[:, k-1] + 37 * f[:, k-2] - 9 * f[:, k-3])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    return(t,y)

# Dibujamos las graficas
figure('Ejercicio 1 apartado a')
tini = perf_counter()
t, y = AB4Sis(ini, fin, fun2, N, y0)
tfin = perf_counter()

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((fin-ini)/N))
print('-----')

subplot(311)
plot(t, y[0], '-*')
title('x frente a t')
subplot(312)
plot(t, y[1], 'r')
title('y frente a t')
subplot(313)
plot(y[0], y[1], 'g')
title('Trayectoria en el plano de fases')
show()

figure('Ejercicio 1 apartado b')
# Dibujemos la frontera de la region de estabilidad absoluta del metodo en el plano complejo
def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

localizar_frontera(array([1, -1, 0, 0, 0]), array([0, 55, -59, 37, - 9])/24) # AB4
title("Fronteras AB")
show()

A = array([[0, 1], [-a, -b]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

rel, iml = -0.5, 3.1225
# Obtenemos h para llevar las raíces a la frontera
h = 0.1
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()

# Como los puntos de los autovalores se encuentran dentro de la región que delimita la 
# frontera de estabilidad, eso indica que los autovalores están dentro de la región 
# de estabilidad absoluta del método. Esto es una buena señal, ya que significa que 
# el método numérico utilizado es estable para los autovalores correspondientes a la
# matriz del sistema.


# Ejercicio 2
# Repita el ejercicio anterior con el método de Adams-Moulton de 4 pasos (AM4):
# yk+1 = yk + h/720 (251 fk+1 + 646 fk − 264k−1 + 106 fk−2 − 19 fk−3), k = 3, 4, . . .
# Use nuevamente el método RK4 para arrancar. Para calcular yk+1 use el método de punto fijo
# tomando como semilla el valor que daría AB4, es decir 
# z0 = yk + h/24 (55 fk − 59 fk−1 + 37 fk−2 − 9 fk−3).
# Detenga la iteración de punto fijo cuando dos términos de la sucesión que se genera disten menos
# de 10−12 o cuando el número de iteraciones llegue a 200. En este último caso, se debe advertir en
# pantalla de que el método de punto fijo tiene problemas para converger.


# Método AM4
# def AM4Sis(a, b, fun, N, y0):
#     y = zeros([len(y0), N+1])
#     t = zeros(N+1)
#     f = zeros([len(y0), N+1])
#     t[0] = a
#     h = (b-a) / N
#     y[:, 0] = y0
#     f[:, 0] = fun(a, y[:, 0])
#     iteraciones = 0


#     ## Metodo de Runge-Kutta de orden 4 Sistemas
#     for k in range(3):
#         k1 = fun(t[k], y[:, k])
#         k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
#         k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
#         k4 = fun(t[k]+h, y[:, k]+h*k3)
#         t[k+1] = t[k]+h
#         y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
#     for k in range(3, N):
#         z = y[k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
#         contador = 0
#         distancia = 1 + 1e-12
#         C = y[k] + h/720 * (251*f[:, k] + 646*f[:, k-1] - 264*f[:, k-2] + 106*f[:, k-3] - 19*f[:, k-4])
#         while distancia > 1e-12 and contador < 100:
#             z_nuevo = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*f[:, k] - 264*f[:, k-1] + 106*f[:, k-2] - 19*f[:, k-3])
#             distancia = norm(z_nuevo - z)
#             z = z_nuevo
#             contador += 1
#         if contador == 100:
#             print("No se ha alcanzado la convergencia en el paso", k)
#         iteraciones = max(iteraciones, contador)
#         y[:, k+1] = z
#         t[k+1] = t[k] + h
#         f[:, k+1] = fun(t[k+1], y[:, k+1])
    
#     return (t, y, iteraciones)


def AM4Sis(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    iteraciones = 0

    ## Método de Runge-Kutta de orden 4 para sistemas
    for k in range(3):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k] + h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
    for k in range(3, N):
        z = y[:, k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
        contador = 0
        distancia = 1 + 1e-12
        C = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*f[:, k] - 264*f[:, k-1] + 106*f[:, k-2] - 19*f[:, k-3])
        while distancia > 1e-12 and contador < 200:
            z_nuevo = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*fun(t[k], y[:, k]) - 264*fun(t[k-1], y[:, k-1]) + 106*fun(t[k-2], y[:, k-2]) - 19*fun(t[k-3], y[:, k-3]))
            distancia = linalg.norm(z_nuevo - z)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("No se ha alcanzado la convergencia en el paso", k)
        iteraciones = max(iteraciones, contador)
        y[:, k+1] = z
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    return (t, y, iteraciones)

# Dibujamos las soluciones
figure("Ejercicio 2a")
tini = perf_counter()
t, y, maxiter = AM4Sis(ini, fin, fun2, N, y0)
tfin = perf_counter()

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((fin-ini)/N))
print('Iteraciones máximas: ' + str(maxiter))
print('-----')

subplot(311)
plot(t, y[0], '-*')
title('x frente a t')
subplot(312)
plot(t, y[1], 'r')
title('y frente a t')
subplot(313)
plot(y[0], y[1], 'g')
title('Trayectoria en el plano de fases')
show()

figure("Ejercicio 2b")
localizar_frontera(array([1, -1, 0, 0, 0]), array([251, 646, -264, 106, -19])/720) # AM4

title("Frontera AB4")
# show()

A = array([[0, 1], [-a, -b]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

rel, iml = -0.5, 3.1225
# Obtenemos h para llevar las raíces a la frontera
h = 0.1
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()



