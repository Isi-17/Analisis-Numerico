from numpy import *
from matplotlib.pyplot import *
from time import perf_counter
from numpy.linalg import eig, norm

# Se considera el problema de Cauchy:
# x'' + 100x' + 25000x = 100*sin(4*pi*t) 
# x(0) = 1,
# x''(0) = 0.

# Reescribe la ecuación como un sistema de ecuaciones de primer orden:

# y[0] = x
# y[1] = x'
def fun(t, y):
    f1 = y[1]   # x' = f1
    f2 = 100*sin(4*pi*t) - 100*y[1] - 25000*y[0] # x'' = f2
    return array([f1, f2])

# Se considera el método de Heun con tablero de Butcher:
# c = array([0, 1])
# A = array([[0, 0], [1, 0]])
# b = array([1/2, 1/2])

# a) Resolver en el intervalo [0, 1] con N = {93, 100, 500} subintervalos. 
# Usando el comando subplot, dibujar en una única figura pero tres gráficas separadas:

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

a = 0
b = 1
y0 = array([1, 0])

N = [93, 100, 500]

figure(1)
for i in range(len(N)):
    t, y = heun_sistema(a, b, fun, N[i], y0)
    subplot(3, 1, i+1)
    plot(t, y[0, :], '-*')
    title(f'N = {N[i]}')
    ylabel('x(t)')
    grid(True)

xlabel('t')
show()

# b) Dibuje la frontera de la region de estabilidad Da del método y compruebe gráficamente si se verifica
# h*lambda_i pertenece a Da  i = 1, 2. Siendo lambda_i los autobalores de la matriz
# Dzf y h = 1/N con N = {93, 100, 500}. 
# Comente si los resultados obtenidos explican algunos aspectos de las gráficas del apartado anterior.

# Dzf = [[0, 1], [-25000, -100]]

Dzf = array([[0, 1], [-25000, -100]])
print("Matriz Dzf:\n", Dzf)
print("Autovalores de Dzf:", eig(Dzf)[0])

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

figure('Ejercicio 1 apartado b frontera de estabilidad')
localizar_frontera(array([1, 1/2, 0]), array([0, 1/2, 0])) # Heun
# Comprobamos si se verifica h*lambda_i pertenece a Da  i = 1, 2

h = [1/93, 1/100, 1/500]
lambda1 = eig(Dzf)[0][0]
lambda2 = eig(Dzf)[0][1]

plot(h[0]*lambda1, h[0]*lambda2, ".", h[0]*lambda1, -h[0]*lambda2, ".")
plot(h[1]*lambda1, h[1]*lambda2, ".", h[1]*lambda1, -h[1]*lambda2, ".")
plot(h[2]*lambda1, h[2]*lambda2, ".", h[2]*lambda1, -h[2]*lambda2, ".")

title("Fronteras Heun")
show()

# Como los puntos de los autovalores se encuentran dentro de la región que delimita la 
# frontera de estabilidad, eso indica que los autovalores están dentro de la región 
# de estabilidad absoluta del método. Esto es una buena señal, ya que significa que 
# el método numérico utilizado es estable para los autovalores correspondientes a la
# matriz del sistema.

# Ejercicio 2
# Repita los dos apartados del ejercicio anterior usando el metodo Adams-Moulton de dos pasos
# z_(k+1) = z_k + h/12 * (5*f(t_k+1, z_k+1) + 8*f(t_k, z_k) - f(t_k-1, z_k-1))
# con N = {66, 100, 500}. Use
# El metodo de Euler para arrancar
# El metodo del punto fijo para calcular z_(k+1), tomando como semilla la aproximacion dada por 
# el metodo de Adams-Bashforth de dos pasos.
# z_(k+1) = z_k + h/2 (3*f(t_k+1, z_k+1) - f(t_k, z_k))
# Detenga las iteraciones del método de punto fijo cuando dos aproximaciones sucesivas verifiquen
# ||z_(k+1) - z_k|| < tol = 1.e-12 o cuando se alcance el numero maximo de iteraciones stop = 500.
# Ponga por pantalla el máximo numero de iteraciones del método de punto fijo que han sido necesarias en cada ejecucion.

# a) Resolver en el intervalo [0, 1] con N = {66, 100, 500} subintervalos.
# Usando el comando subplot, dibujar en una única figura pero tres gráficas separadas:

# Vamos a definir el metodo AM2 para este problema
def AM2pfSis(a,b,fun,N,y0):
        
        tol = 1.e-12
        stop = 500
        
        y = zeros([len(y0),N+1])
        t = zeros(N+1)
        f = zeros([len(y0),N+1])
        t[0] = a
        h = (b-a)/float(N) 
        y[:,0] = y0
        f[:,0] = fun(a,y0)
        maxiter = 0
        for k in range(1):  ### Arranco con Euler
            t[k+1] = t[k]+h
            y[:, k+1] = y[:, k]+h*fun(t[k], y[:, k])
            f[:,k+1] = fun(t[k+1],y[:,k+1])
        for k in range(1,N):
            Ck = y[:,k] +h/2*(3*f[:,k] - f[:,k-1])
            t[k+1] = t[k]+h
            cont = 0
            z = y[:,k]
            error = tol +1
            while(error >=tol and cont < stop):
                znew = h*fun(t[k+1], z) +Ck
                error = max(abs(znew - z))
                z = znew
                cont += 1
            maxiter = max(maxiter,cont)
            if (cont == stop):
                print('El método de punto fijo no ha convergido tras '+str(stop)+' iteraciones para N = '+str(N)+'.')
            y[:,k+1] = z
            f[:,k+1] = fun(t[k+1],y[:,k+1])
        return(t,y,maxiter)


a = 0
b = 1
y0 = array([1, 0])
N = [66, 100, 500]
figure('Ejercicio 2 apartado a')
for i in range(len(N)):
    t, y, maxiter = AM2pfSis(a, b, fun, N[i], y0)
    subplot(3, 1, i+1)
    plot(t, y[0, :], '-*')
    title(f'N = {N[i]}')
    ylabel('x(t)')
    grid(True)
    if i == len(N)-1:
        xlabel('t')
    print(f'El máximo número de iteraciones del método de punto fijo para N = {N[i]} es {maxiter}.')
show()

# 2.b) Dibuje la frontera de la region de estabilidad Da del método y compruebe gráficamente si se verifica
# h*lambda_i pertenece a Da  i = 1, 2. Siendo lambda_i los autobalores de la matriz
# Dzf y h = 1/N con N = {66, 100, 500}. 
# Comente si los resultados obtenidos explican algunos aspectos de las gráficas del apartado anterior.
# Dzf = [[0, 1], [-25000, -100]]

# Dibujemos la frontera de la region de estabilidad absoluta del metodo en el plano complejo
figure('Ejercicio 2 apartado b frontera de estabilidad')
localizar_frontera(array([1, -1, 0]), array([5, 8, -1])/12) # AM2
# Comprobamos si se verifica h*lambda_i pertenece a Da  i = 1, 2

h = [1/66, 1/100, 1/500]
lambda1 = eig(Dzf)[0][0]
lambda2 = eig(Dzf)[0][1]

plot(h[0]*lambda1, h[0]*lambda2, ".", h[0]*lambda1, -h[0]*lambda2, ".")
plot(h[1]*lambda1, h[1]*lambda2, ".", h[1]*lambda1, -h[1]*lambda2, ".")
plot(h[2]*lambda1, h[2]*lambda2, ".", h[2]*lambda1, -h[2]*lambda2, ".")

title("Fronteras Heun")
show()

# Como los puntos de los autovalores se encuentran dentro de la región que delimita la 
# frontera de estabilidad, eso indica que los autovalores están dentro de la región 
# de estabilidad absoluta del método. Esto es una buena señal, ya que significa que 
# el método numérico utilizado es estable para los autovalores correspondientes a la
# matriz del sistema.




