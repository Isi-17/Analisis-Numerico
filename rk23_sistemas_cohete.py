from pylab import *
from time import perf_counter

g = 9.81
M = 7.5
alpha = 0.02
mfuel0 = 7.5


    
def fun(t,z):
    #z[0]:= x; z[1]:=y; z[2]:=velocidad, z[3]:=angulo, z[4]:=masa de combustible 
    
    T = T0*(z[2] >0 )
    m = M + z[2]
    
    f0 = z[1]
    f1 = -g + T/m - C/m*z[1]*abs(z[1]) + alpha*T/m*z[1]
    f2 = -alpha*T
    
    return array([f0,f1,f2])


def rk23_sis_rocket(a, z0, h0, tol):
    """ Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0 """
    

    
    maxstep = 10000 #el programa se detiene si se han dado 10000 pasos para evitar bucles infinitos

    hmin = 1.e-5 # paso de malla minimo
    hmax = .1 # paso de malla maximo
    
    
    # Coeficientes RK
    q = 3 # etapas del metodo
    p = 2 # orden del método menos preciso
    
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
    
    # Inicializacion de variables
    t = array ([a]) # nodos
    z = zeros([len(z0),1])
#    z = z0.reshape(len(z0), 1) # alterativa para poner el vector de c.i.como columna
    z[:,0] = z0
    h = array ([h0]) # pasos de malla
    K = zeros ([len(z0), q])
    k = 0 # contador de iteraciones
    nonstop = True
    
    while (nonstop and k < maxstep):
        for i in range(q):
            K[:,i] = fun(t[k]+C[i]*h[k], z[:,k]+h[k]*dot(A[i,:], transpose(K)))
        incrlow = dot(B, transpose(K))
        incrhigh = dot(BB, transpose(K))
        error = h[k]*(incrhigh-incrlow) # estimacion del error
        z = column_stack((z, z[:,k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])
        hnew = 0.9*h[k]*abs(tol/norm(error, inf))**(1./(p+1))
        hnew = min(max(hnew, hmin), hmax)
        h = append(h, hnew)
        k += 1
        if z[0,k] <= 0:
            nonstop = False
    if k == maxstep:
        print('maximo numero de pasos alcanzados')
    return (t, z, h)

a=0.
h0 = 0.05

z0 = 0.0
v0 = 50.0
tol = 1.e-4


z0=array([z0,v0,mfuel0])

print('-----')
print('1. Sin motor ni rozamiento')
T0 = 0
C = 0
tini = perf_counter()
t,z,h=rk23_sis_rocket(a, z0, h0, tol)
tfin = perf_counter()

hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla


print("No. de iteraciones: " + str(len(t)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print("Tiempo de caida: "+str(sum(h[:-1])))
print('Máxima altura: ' + str(max(z[0,:])))

subplot(211)
plot(t,z[0])
xlabel('t')
ylabel('z')
title('Altura del cohete')
subplot(212)
plot(t[:-1],h[:-1])
title('Pasos usados')
print('-----')

print('2. Sin motor con rozamiento')
T0 = 0
C = 0.02
tini = perf_counter()
t,z,h=rk23_sis_rocket(a, z0, h0, tol)
tfin = perf_counter()

hn = min(h[:-1]) # minimo valor utilizado del paso de malla
hm = max(h[:-1]) # maximo valor utilizado del paso de malla


print("No. de iteraciones: " + str(len(t)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print("Tiempo de caida: "+str(sum(h[:-1])))
print('Máxima altura: ' + str(max(z[0,:])))

subplot(211)
plot(t,z[0])
subplot(212)
plot(t[:-1],h[:-1])

print('-----')


print('3. Con motor con rozamiento')
T0 = 50
C = 0.02
tini = perf_counter()
t,z,h=rk23_sis_rocket(a, z0, h0, tol)
tfin = perf_counter()

hn = min(h[:-1]) # minimo valor utilizado del paso de malla
hm = max(h[:-1]) # maximo valor utilizado del paso de malla

print("No. de iteraciones: " + str(len(t)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print("Tiempo de caida: "+str(sum(h[:-1])))
print('Máxima altura: ' + str(max(z[0,:])))

subplot(211)
plot(t,z[0])
xlabel('t')
ylabel('z')
subplot(212)
plot(t[:-1],h[:-1])
print('-----')

subplot(211)
legend(['C=0, T=0', 'C= 0.2, T = 0', 'C = 0.2, T= 50'])
subplot(212)
legend(['C=0, T=0', 'C= 0.2, T = 0', 'C = 0.2, T= 50'])

