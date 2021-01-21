%% Experimentos para cálculo de parámetros
clc;
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
m2 = 3.19e-3; % Kg
tau_2 = 2e-3; % Nm
t_2 = 0.6592; % rad
% Se pone escalon de tau_2 = g*m2*l2*sin(t_2)
l2 = tau_2/(g*m2*sin(t_2));
%% C) J2 - Momento de inercia 2
% Se usa el expermiento C 
p = 721.604e-3;
J2 = l2*m2*g*((p/(2*pi))^2);

%% Linealizacion de ecuaciones de estado
%b1 y b2: damping coefficient
b1 = 1e-6; %cambiar
b2 = 0.001e-3; %cambiar
L1 = 103.5e-3;  % Abrimos las dos piezas con un visualizador
[theta,A,B] = linealizacion(m1,m2,l1,l2,L1,J1,J2,tau_1,b1,b2,g);
C = [0 1 0 0];
%% Realimentacion de estados (No estaria andando, se queda oscilando)
% Pruebo sin accion integral
% pole = -10;
% K = 180*acker(A,B,[pole-1;pole-2;pole-3;pole-4])/pi;
% Con accion integral
% Ba = [B;0];
% Aa = [A zeros(4,1);-C 0];
% Kt = acker(Aa,Ba,[-10;-10;-10;-10;-10]);
% K = Kt(1:end-1)
% Ki = -Kt(end)
% Prueba con LQR (Ahora anda)
Q = diag([10 10 10 10]);
R = 5;
[K,S,e] = lqr(A,B,Q,R);
% Prueba con LQI
% Aa = [A zeros(4,1);-C 0];
% Ba = [B;0];
Qi = diag([2 2 10 10 10]);
R = 7;
sys = ss(A,B,C,0);
[Kt,St,et] = lqi(sys,Qi,R);
Kt
%% LoopShaping
% Primero sacamos la transferencia SIMO
% C = [1 0 0 0];
% [n,d] = ss2tf(A,B,C,0);
% P1 = tf(n,d);
% C = [0 1 0 0];
% [n,d] = ss2tf(A,B,C,0);
% P2 = tf(n,d);
% % Disenio controlador C2
% s = tf('s');
% C2 = (1/s);