close all; clear all; clc;

%% Setting up constants

% Constants
params.mq = 0.5;
params.mp = 1;
params.L = 0.5;
params.e3 = [0;0;1];
params.g = -9.81;

% Controller & aero Constants
params.pos_natFreq = 2.0;                 % [rad/s] natural frequency for position loop
params.pos_damping = 0.7;                 % damping ratio
params.maxProperAcc = 13;                 % [m/s^2] max proper acceleration
params.minVerticalProperAcceleration = -1; % [m/s^2] min vertical proper acceleration
params.linDragCoeffB = [0.1286; 0.1286; 0.1286]; % body-frame linear drag [N/(m/s)]
params.payloadLinDragCoeffB = [0; 0; 0];    % payload linear drag [N/(m/s)] in world frame

%% Analysis

mq = params.mq;
mp = params.mp;
e3 = params.e3;
g = params.g;
L = params.L;

% Quadcopter drag in world frame:
v_b = R_wb.' * dq;
K = diag(params.linDragCoeffB);   % [N/(m/s)]
F_drag_b = -K * v_b;
F_drag_w = R_wb * F_drag_b;

% Payload drag in world frame
Kp_payload = diag(params.payloadLinDragCoeffB);
F_drag_p = -Kp_payload * dp;

% Initial conditions
rhat = [0;0;-1];

% s_initial = [q0;dq0];
% s0 = s_initial;
% 
% sp_initial = [qp0; dqp0];
% sp0 = sp_initial;

A = [mq * eye(3), zeros([3,3]), -rhat;
     zeros([3,3]), mp * eye(3), +rhat;
     -rhat.', rhat.', 0];

B = [ f + F_drag_w + e3 * mq * g;
      e3 * mp * g + F_drag_p;
      -norm(dp - dq)^2 / L ];

X = A\B;