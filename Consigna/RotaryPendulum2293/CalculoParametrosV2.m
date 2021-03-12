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
J1 = 2*pi*tau_1/(tdd_1);
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
B = B/1000; % Compensar el N*mm

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
Lobs = place(Aod',Cod',[poles2d(1) poles2d(1)].^20);
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
[n,d] = ss2tf(A,B,[0 1 0 0],0);
Ga = minreal(zpk(tf(n,d)),0.01); % Sale con beta

% Para beta armamos un PI
C2 = minreal((s+9.518)/(1+(s/50)),0.01);
Lazo2 = minreal(C2*Ga,0.01);
figure();
bode(Lazo2);
Kl2 = 1/db2mag(-12);
Lazo2 = Lazo2*Kl2;
figure();
nyqlog(Lazo2);

C2_aux = C2;
C2 = -C2*Kl2;
C2 = ss(C2);

C2.u='beta';
C2.y='ut';

% Interconectamos C2 con el Sistema Completo:
% Nos queda el sistema con C2 enganchado que realimenta el angulo "beta"

Gss=ss(A,B,[1 0 0 0;0 1 0 0],[0;0]);
Gss.u='u';
Gss.y='y';
Sum = sumblk('u = ut + up');
Sys=ss([],[],[],[1 0]);
Sys.u='y';
Sys.y='alfa';
Sys2=ss([],[],[],[0 1]);
Sys2.u='y';
Sys2.y='beta';

Gb=connect(Gss,C2,Sum,Sys,Sys2,'up','alfa');
Gb = zpk(Gb);
Gb = minreal(Gb,0.01);

% Para alfa armamos un PID
close all;
C1p = minreal((1/s)*(s+9.518)*(s^2 + 40.52*s + 631)*(s+8.696)/(6.94*(s+50)*(s+8.719)),0.01);
zpk(minreal(C1p*zpk(Gb),0.01));
C1pp = -((s+0.15)^2/(1+(s/40))^3);
C1 = minreal(C1p*C1pp);
Lazo1 = minreal(C1*Gb);
figure();
bode(Lazo1);
Kl1 = 1/db2mag(40);
Lazo1 = Lazo1*Kl1;
figure();
nyqlog(Lazo1);
close all;

% Tiempo Discreto
% Polo mas lejos: -50 => 8Hz
% Para cumplir Nyquist al menos muestrear a 16Hz => 6.25ms
Ts2 = 1e-3; % Funciona bien, no oscila
C2d = c2d(C2_aux,Ts2,'zoh');
C1d = c2d(C1,Ts2,'zoh');