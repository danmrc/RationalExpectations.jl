#Tests of Klein and Sims Algorithms for solving Rational Expectations Model
#
#Based on Chapter 3 of Gali

using RationalExpectations
using Plots

#Parameters

sigma = 1
phi_pi = 1.5
phi_y = 0.5/4
beta = 0.99
phi = 1
alph = 1/3
ep = 6
theta = 2/3
rho_v = 0.5

Theta = (1-alph)/(1-alph+alph*ep)
lambda = (1-theta)*((1-beta*theta)/theta)*Theta
kappa = lambda*(sigma+(phi+alph)/(1-alph))

#True solution

Lambda_v = 1/((1-beta*rho_v)*(sigma*(1-rho_v)+phi_y)+kappa*(phi_pi - rho_v))
true_y(v) = -(1-beta*rho_v)*Lambda_v*v
true_pi(v) = -kappa*Lambda_v*v

true_path = zeros(21,3)
initial_shock = 0.25
shock = zeros(21)
shock[2] = initial_shock

for j = 2:21
    true_path[j,1] = rho_v*true_path[j-1] + shock[j]
    true_path[j,2] = true_y(true_path[j,1])
    true_path[j,3] = true_pi(true_path[j,1])
end

plot(true_path[2:21,1])
plot(true_path[2:21,2])
plot(4*true_path[2:21,3])

#Klein

A = [[1 0 0 0];[0 1 0 0]; [0 0 beta 0]; [0 -1 1 sigma]]
B = [[rho_v 0 0 0];[1 0 phi_pi phi_y];[0 0 1 -kappa];[0 0 0 sigma]]

C = [1;0; 0; 0]

k0 = [0.25;0]

t=12

choque = [0.; 0;0; 0]

x = klein(A,B,C,t,k0,choque,[3 4])
plot(x[:,1])
plot(4*x[:,2])
plot(4*x[:,3])
plot(x[:,4])

#Sims

G0 = [[1 0 0 0];[0 1 0 0];[0 -1 sigma 1];[0 0 0 beta]]
G1 = [[rho_v 0 0 0];[1 0 phi_y phi_pi];[0 0 sigma 0];[0 0 -kappa 1]]
Psi = [1 0 0 0]'
Pi = [[0 phi_pi 0 1];[0 phi_y sigma -kappa]]'

matrizes = sims(G0,G1,Pi,Psi)

resul = irf(real(matrizes.Theta1),real(matrizes.Theta2),20,0.25)

plot(0:20,resul[:,1])
plot(4*resul[:,2])
plot(resul[:,3])
plot(4*resul[:,4])


#Indeterminate case

gamma = 0.8
kappa = 0.5
beta = 0.9
rho_t = 0.8

G0 = [[1 0 0 0 0];[-1/gamma 1 -1/gamma 0 0];[0 0 beta 0 0];[0 0 0 1 0];[0 0 0 0 1]]
G1 = [[rho_t 0 0 0 0];[0 1 0 0 0];[0 -kappa 1 0 0];[0 1 0 0 0];[0 0 1 0 0]]
PI = [[0 0];[1 0];[-kappa 1];[1 0];[0 1]]
PSI = [[1 0 ];[0 -1];zeros(3,2)]

mm = sims(G0,G1,PI,PSI)

irf(mm.theta1,mm.theta2,20)

ttest = zeros(12,4)

for j = 1:12
    ttest[j,:] = (C + (matrizes.theta1^(j-1))*matrizes.theta2)*0.25
end
