function [theta,A,B] = linealizacion(m1,m2,l1,l2,L1,J1,J2,tau_1,b1,b2,g)
%% Ecuaciones de estado linealizadas 
% Se utiliza la ecuacion (36) del paper "On the Dynamics of the Furuta Pendulum"
% j0s = J1 + m1*(l1^2) + m2*(L1^2);
% j2s = J2 + m2*(l2^2);

j1s = J1;
j0s = j1s + m2*(L1^2);
j2s = J2;

A31 = 0;
A32 = (g*L1*(m2*l2)^2) / (j0s*j2s - (m2*L1*l2)^2);
A33 = (-b1*j2s) / (j0s*j2s - (m2*L1*l2)^2); 
A34 = (-b2*m2*l2*L1) / (j0s*j2s - (m2*L1*l2)^2);
A41 = 0;
A42 = (g*m2*l2*j0s) / (j0s*j2s - (m2*L1*l2)^2);
A43 = (-b1*m2*l2*L1) / (j0s*j2s - (m2*L1*l2)^2);
A44 = (-b2*j0s) / (j0s*j2s - (m2*L1*l2)^2);
B31 = j2s / (j0s*j2s - (m2*L1*l2)^2);
B41 = (m2*L1*l2) / (j0s*j2s - (m2*L1*l2)^2);


syms t_1 t_2 td_1 td_2;
tdd_1 = [A31 A32 A33 A34]*[t_1; t_2; td_1; td_2] + B31*tau_1;
tdd_2 = [A41 A42 A43 A44]*[t_1; t_2; td_1; td_2] + B41*tau_1;

theta = [td_1; td_2; tdd_1; tdd_2];
A = [0 0 1 0;0 0 0 1;A31 A32 A33 A34;A41 A42 A43 A44];
B = [0;0;B31;B41];
end








