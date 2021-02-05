%% Experimentos para cálculo de parámetros
clear all;
close all;
clc;
run("RotaryPendulum_PARAM.m");
g = 9.81;
m1 = 0.01835; % Suma de las masas
m2 = 0.00575; % Suma de las masas
%% B) J1 - Momento de inercia 1
tau_1 = 1;
% Del grafico sale la aceleracion
tdd_1 = 5.71e4;
% Calculo del J1
J1 = tau_1/(tdd_1);
%% D) l2 - Distancia centro de masa 2
% Haciendo tau_1 = 0
tau_2 = 2e-3; % Nm
t_2 = 0.6592; % rad
% Se pone escalon de tau_2 = g*m2*l2*sin(t_2) L2 = 94mm + 79mm = 173mm
l2 = tau_2/(g*m2*sin(t_2));
%% C) J2 - Momento de inercia 2
p = 721.604e-3;
% Sabiendo la oscilacion
J2 = l2*m2*g*((p/(2*pi))^2);

%% Linealizacion de ecuaciones de estado
%b1 y b2: damping coefficient
b1 = 1e-6; %cambiar
b2 = 0.001e-3; %cambiar
L1 = 103.5e-3;  % Abrimos las dos piezas con un visualizador
                % UltimakerCura
[theta,A,B] = linealizacion(m1,m2,0,l2,L1,J1,J2,tau_1,b1,b2,g);
C = [1 0 0 0];
B = B/800;

%% Realimentacion de estados
% Tiempo Continuo
% 1) LQR sin accion integral
Q1 = diag([1 7 1 7]);
R1 = 1;
[K1,S1,e1] = lqr(A,B,Q1,R1);
poles1 = eig(A-B*K1);
% 2) LQI
sys = ss(A,B,C,0);
Q2 = diag([1 7 1 7 0.05]);
R2 = 1;
[K2,S2,e2] = lqi(sys,Q2,R2);
Aa = [A zeros(4,1);-C 0];
Ba = [B;0];
poles2 = eig(Aa-Ba*K2);
% 3) LQI con observador
% Hacemos modo reducido
% - Mido alfa y beta
% - Estimo alfa_d y beta_d
Aw = A;
Bw = B;
Aaa = Aw(1:2,1:2);
Aab = Aw(1:2,3:4);
Aba = Aw(3:4,1:2);
Abb = Aw(3:4,3:4);
Ba = Bw(1:2);
Bb = Bw(3:4);
% Para el calculo de L
Ao = Abb;
Co = Aab;
Lobs = place(Ao',Co',10*[poles2(1) poles2(1)]);
Ke = Lobs';
poles3_obs = eig(Ao-Ke*Co);
% Matrices equivalentes para simulink
A_h = Abb - Ke*Aab;
B_h = A_h*Ke + Aba - Ke*Aaa;
F_h = Bb - Ke*Ba;
C_h = [0 0;0 0;1 0;0 1];
D_h = [1 0;0 1;Ke];

% Tiempo Discreto
Ts = 10e-3;
sys = ss(A,B,C,0);
sysD = c2d(sys,Ts,'zoh');
Ad = sysD.a;
Bd = sysD.b;
Cd = sysD.c;
% 1d) LQR sin accion integral
Q1d = diag([1 7 1 7]);
R1d = 1;
[K1d,S1d,e1d] = dlqr(Ad,Bd,Q1d,R1d);
poles1d = eig(Ad-Bd*K1d);
% 2d) LQI
Q2d = diag([1 7 1 7 0.5]);
R2d = 1;
[K2d,S2d,e2d] = lqi(sysD,Q2d,R2d);
Aad = [Ad zeros(4,1);-Cd*Ts 1];
Bad = [Bd;0];
poles2d = eig(Aad-Bad*K2d);
% 3d) LQI con observador
% Hacemos modo reducido
% - Mido alfa y beta
% - Estimo alfa_d y beta_d
Aw = Ad;
Bw = Bd;
Aaa = Aw(1:2,1:2);
Aab = Aw(1:2,3:4);
Aba = Aw(3:4,1:2);
Abb = Aw(3:4,3:4);
Ba = Bw(1:2);
Bb = Bw(3:4);
% Para el calculo de L
Aod = Abb;
Cod = Aab;
Lobs = place(Aod',Cod',[poles2d(1) poles2d(1)].^10);
Ked = Lobs';
poles3d_obs = eig(Aod-Ked*Cod);
% Matrices equivalentes para simulink
Ad_h = Abb - Ked*Aab;
Bd_h = Ad_h*Ked + Aba - Ked*Aaa;
Fd_h = Bb - Ked*Ba;
Cd_h = [0 0;0 0;1 0;0 1];
Dd_h = [1 0;0 1;Ked];

%% Loop Shaping
s = tf('s');
[n,d] = ss2tf(A,B,[1 0 0 0],0);
P1 = minreal(zpk(tf(n,d)),0.001);
[n,d] = ss2tf(A,B,[0 1 0 0],0);
Ga = minreal(zpk(tf(n,d)),0.001); % Sale con beta
pade = (1-(s*Ts/4))/(1+(s*Ts/4));
Gb = minreal(P1/Ga,0.001);
% Tal que Ga*Gb = P1
% Para beta armamos un PI
Ga_mp = minreal(Ga*(s-10.77)/(s+10.77),0.001);
C2 = minreal((1/Ga_mp)/((1+(s/50))^2),0.001);
Lazo2 = minreal(C2*Ga,0.001);
figure();
bode(Lazo2);
Kl2 = 1/db2mag(-35);
Lazo2 = Lazo2*Kl2;
figure();
nyqlog(Lazo2);
% Para alfa armamos un PID
close all;
T2 = minreal(1-(1/(1+Lazo2)),0.001);
Gbeq = minreal(Gb*T2,0.001);
Gbeq_mp = minreal(Gbeq*(s+8.696)/((s-8.696)),0.001);
C1 = minreal((-1)*(1/Gbeq_mp)/((1+(s/0.5))^2*(1+(s/1))*(1+(s/50))),0.001);
Lazo1 = minreal(C1*Gbeq,0.001);
figure();
bode(Lazo1);
Kl1 = 1/db2mag(-10);
Lazo1 = Lazo1*Kl1;
figure();
nyqlog(Lazo1);
close all;
T1 = minreal(1-(1/(1+Lazo1)),0.001);
eig(T1)