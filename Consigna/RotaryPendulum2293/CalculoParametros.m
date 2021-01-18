%% Experimentos para cálculo de parámetros
g = 9.81;
%% B) J1 - Momento de inercia 1
tau_1 = 1;
% Del grafico sale la aceleracion
tdd_1 = 5.71e4;
% Calculo del J1
J1 = tau_1/tdd_1;
%% C) l1 - Distancia centro de masa 1
p = 721.604e-3;
m1 = 0.0198; % Suma de las masas
% Sabiendo la oscilacion
l1 = J1/(m1*g*((p/(2*pi))^2));
%% D) l2 - Distancia centro de masa 2
% Haciendo tau_1 = 0
m2 = 3.5e-3; % Kg
tau_2 = 2e-3; % Nm
t_2 = 0.6592; % rad
% Se pone escalon de tau_2 = g*m2*l2*sin(t_2)
l2 = tau_2/(g*m2*sin(t_2));
%% D) J2 - Momento de inercia 2
% Se usa el expermiento D sin dumping y sin gravedad
tdd_2 = 46.6;
J2 = tau_2/tdd_2;