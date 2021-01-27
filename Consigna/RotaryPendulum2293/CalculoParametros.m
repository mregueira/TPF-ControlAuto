%% Experimentos para cálculo de parámetros
clear all;
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
%% Realimentacion de estados
C = [1 0 0 0];
% Pruebo sin accion integral (Funciona OK, todo +)
pole1 = -2.8+3i;
pole2 = -2.8-3i;
Kn = acker(A,B,130*[pole1;pole2;-0.03;-0.02]);
% Con accion integral (Funciona OK, todo +)
Ba = [B;0];
Aa = [A zeros(4,1);-C 0];
Kt = acker(Aa,Ba,5*[pole1;pole2;-0.03;-0.02;-3000]);
K = Kt(1:end-1)
Ki = -Kt(end)
% Prueba con LQR (funciona OK)
% Q = diag([10 10 10 10]);
% R = 5;
% [K,S,e] = lqr(A,B,Q,R);
% Prueba con LQI (funciona OK)
% Qi = diag([2 2 10 10 10]);
% R = 7;
% sys = ss(A,B,C,0);
% [Kt,St,et] = lqi(sys,Qi,R);
% Kt;

% En discreto (Funciona OK, todo +, 1ms)
% - Uso las ganancias continuas
% - Introducir ZOH
% - Cambiar integrador por Forward
Ts = 1e-3;
% sysD = c2d(sys,Ts,'zoh');
% Aa = [sysD.a zeros(4,1);-Ts*C 1];
% Ba = [sysD.b;0];
% mul = 8;
% p1 = exp(pole1*Ts*mul);
% p2 = exp(pole2*Ts*mul);
% p3 = exp(-0.08*Ts*mul);
% p4 = exp(-0.07*Ts*mul);
% p5 = exp(-5000*Ts*mul);
% Kd = acker(Aa,Ba,[p1 p2 p3 p4 p5]);
% Kd_f = Kd(1:end-1)
% Kd_i = -Kd(end)
Lobs = acker(A',C',130*[pole1;pole2;-0.03;-0.02])
% Qi = diag([1 1 20 10 7]);
% R = 0.7;
% [Kd,Sd,ed] = dlqr(AaD,BaD,Qi,R);
% Kd
%% LoopShaping
% % Primero sacamos la transferencia
% C = [1 0 0 0];
% [n,d] = ss2tf(A,B,C,0);
% P1 = tf(n,d);
% P_poles = pole(P1)
% P1_zeros = zero(P1)
% C = [0 1 0 0];
% [n,d] = ss2tf(A,B,C,0);
% P2 = tf(n,d);
% P2_zeros = zero(P2);
% % Simplificado
% s = tf('s');
% P1_s = 2.677e4*(s+6.4)*(s-6.4)/(s*(s+7.7)*(s-7.7)*(s+0.02));
% P2_s = 1.185e4*s/((s+7.7)*(s-7.7)*(s+0.02));
% % Disenio controlador C2
% C2 = (s+7.7)^2 * (s+0.02) /(s * (1+(s/100))^2);
% L2 = P2_s*C2;
% figure;
% bode(L2)
% k2 = 1/db2mag(75);
% C2 = C2*k2;
% L2 = P2_s*C2;
% figure;
% bode(L2)
% figure;
% nyqlog(L2)
% % Disenio controlador C1
% close all;
% clc
% Peq_s = minreal(C2*P1_s/(1+C2*P2_s),0.0001);
% C1 = s^2 * (s+0.03) * (s^2+189.4*s+2.89e4) / ((s+6.4)^2 * (s+7.7)*((s/100)+1)^2);
% L1 = Peq_s*C1;
% figure;
% bode(L1)
% k1 = 1/db2mag(95);
% C1 = C1*k1;
% L1 = Peq_s*C1;
% figure;
% bode(L1)
% figure;
% nyqlog(L1)