clear
clc

syms x th real
syms dx ddx dth ddth real
syms m_c m_l l g real

% m_c = 0.10;
% m_l = 0.01;
% l = 1;
% g = 9.81;

q   = [x; th];      % Generalized Coordinates
dq  = [dx; dth];    % Velocities

T = 1/2*m_c*dx^2 + 1/2*m_l*l^2*dth^2 + 1/2*m_l*dx^2 + m_l*l*dx*dth*cos(th);
U = -m_l*l*g*cos(th);

L = T - U;

q_act = [1; 1];

[D, C, G, B] = LagrangianDynamics(T, U, q, dq, q_act) ;

matlabFunction(D, 'File', 'auto_D')
matlabFunction(C, 'File', 'auto_C')
matlabFunction(G, 'File', 'auto_G')
matlabFunction(B, 'File', 'auto_B')