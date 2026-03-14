close all; clear all; clc;
colororder(groot, 'gem');


g = 9.81;

% Constants
params.mq = 0.7;
params.mp = 0.7;
params.L = 0.5;
params.e3 = [0;0;1];
params.g = -9.81;

q0 = [0;
      0;
      2;];

dq0 = [0;
       0;
       0;];

qp0 = [0;
      0;
      q0(3) - params.L;];

dqp0 = [0;
       0;
       0;];

% Controller & aero Constants
params.pos_natFreq = 2.0;                 % [rad/s] natural frequency for position loop
params.pos_damping = 0.7;                 % damping ratio
params.maxProperAcc = 13;                 % [m/s^2] max proper acceleration
params.minVerticalProperAcceleration = -1; % [m/s^2] min vertical proper acceleration
params.linDragCoeffB = [0.1286; 0.1286; 0.1286]; % body-frame linear drag [N/(m/s)]

payloadLinDrag = 1.225*3.14*(0.067/2)^2*0.6*3;

params.payloadLinDragCoeffB = [payloadLinDrag; payloadLinDrag; payloadLinDrag];    % payload linear drag [N/(m/s)]

% Initial conditions
rhat = [0;0;-1];

s_initial = [q0;dq0];
s0 = s_initial;

sp_initial = [qp0; dqp0];
sp0 = sp_initial;

% T = params.mp * 9.81;

dt = 0.002;
t_steps = 0:dt:10;

s_arr = zeros(6, length(t_steps));
s_arr(:, 1) = s0; % Initialize the first state in the array

sp_arr = zeros(6, length(t_steps));
sp_arr(:, 1) = sp0; % Initialize the first state in the array

T_arr = zeros(1, length(t_steps));
% T_arr(1,1) = T;

R_arr = zeros(3, length(t_steps));
R_arr(:,1) = rhat;

for i = 2:length(t_steps)
    t = t_steps(i);
    [s1, sp1, T1, rhat] = simulate(t, dt, s0, sp0, rhat, params);
    s_arr(:, i) = s1; % Store the state in the array for each time step
    sp_arr(:, i) = sp1; % Store the state in the array for each time step

    T_arr(1, i) = T1;
    T = T1;

    R_arr(:, i) = rhat;

    s0 = s1; % Update the state for the next iteration
    sp0 = sp1;

end

figure(1); hold on;
title("Quadcopter Position over time (mass: " + params.mp + "kg)")
plot(t_steps, s_arr(1,:), 'r', 'DisplayName', 'x', 'LineWidth', 2)
plot(t_steps, s_arr(2,:), 'g', 'DisplayName', 'y', 'LineWidth', 2)
plot(t_steps, s_arr(3,:), 'b', 'DisplayName', 'z', 'LineWidth', 2)
legend;

figure(2); hold on;
title("Payload Position over time (mass: " + params.mp + "kg)")
plot(t_steps, sp_arr(1,:), 'r', 'DisplayName', 'x', 'LineWidth', 2)
plot(t_steps, sp_arr(2,:), 'g', 'DisplayName', 'y', 'LineWidth', 2)
plot(t_steps, sp_arr(3,:), 'b', 'DisplayName', 'z', 'LineWidth', 2)
legend;


figure(3);
subplot(2, 1, 1); hold on;
title("Quadcopter Position over time (mass: " + params.mq + "kg)")
plot(t_steps, s_arr(1,:), 'r', 'DisplayName', 'x', 'LineWidth', 2)
plot(t_steps, s_arr(2,:), 'g', 'DisplayName', 'y', 'LineWidth', 2)
plot(t_steps, s_arr(3,:), 'b', 'DisplayName', 'z', 'LineWidth', 2)
legend;

subplot(2, 1, 2); hold on;
title("Payload Position over time (mass: " + params.mp + "kg)")
plot(t_steps, sp_arr(1,:), 'r', 'DisplayName', 'x', 'LineWidth', 2)
plot(t_steps, sp_arr(2,:), 'g', 'DisplayName', 'y', 'LineWidth', 2)
plot(t_steps, sp_arr(3,:), 'b', 'DisplayName', 'z', 'LineWidth', 2)
legend;

figure(4); hold on;
plot(t_steps(2:end), T_arr(1,2:end), 'r', 'LineWidth', 2)
title('Tension Force over time')

figure(5); hold on;
title('rhat over time')
plot(t_steps, R_arr(1,:), 'r', 'DisplayName', 'x', 'LineWidth', 2)
plot(t_steps, R_arr(2,:), 'g', 'DisplayName', 'y', 'LineWidth', 2)
plot(t_steps, R_arr(3,:), 'b', 'DisplayName', 'z', 'LineWidth', 2)

% animate([s_arr;sp_arr], dt)

% figure(2)
% plot3(s_arr(1,:), s_arr(2,:), s_arr(3,:), '.')

%% Helper functions

function [s1, sp1, T1, new_rhat] = simulate(t, dt, s0, sp0, rhat, params)

mq = params.mq;
mp = params.mp;
e3 = params.e3;
g = params.g;
L = params.L;

q = s0(1:3);
dq = s0(4:6);

p = sp0(1:3);
dp = sp0(4:6);

[f, R_wb] = u_input(params, t, q, dq);

% Body-frame linear drag:
v_b = R_wb.' * dq;
K = diag(params.linDragCoeffB);   % [N/(m/s)]
F_drag_b = -K * v_b;
F_drag_w = R_wb * F_drag_b;

% Payload drag in world frame
Kp_payload = diag(params.payloadLinDragCoeffB);
F_drag_p = -Kp_payload * dp;

A = [mq * eye(3), zeros([3,3]), -rhat;
     zeros([3,3]), mp * eye(3), +rhat;
     -rhat.', rhat.', 0];

B = [ f + F_drag_w + e3 * mq * g;
      e3 * mp * g + F_drag_p;
      -norm(dp - dq)^2 / L ];

X = A\B;

ddq = X(1:3);
ddp = X(4:6);
T1 = X(7);

new_dq = dq + dt*ddq;
new_dp = dp + dt*ddp;

new_q = q + dt*new_dq;
new_p = p + dt*new_dp;

new_rhat = (new_p-new_q)/norm(new_p-new_q);

sp1 = [new_p; new_dp];
s1 = [new_q; new_dq];

end

function [f, R] = u_input(params, t, p, pdot)
    % Position controller that mirrors the C++ QuadcopterController::Run:
    %   - second-order (wn, zeta) position loop -> desired acceleration
    %   - add gravity to get proper acceleration
    %   - saturate proper acceleration magnitude and vertical component
    %   - construct attitude so body-z aligns with thrust direction

    m = params.mq + params.mp;
    g = 9.81;

    p_des = desired_position(t);
    pdot_des = [0;0;0];
    desAcc = [0;0;0];

    wn = params.pos_natFreq;
    zeta = params.pos_damping;

    Kp = wn^2 * eye(3);
    Kd = 2 * zeta * wn * eye(3);

    cmdAcc = Kp*(p_des - p) + Kd*(pdot_des - pdot) + desAcc;
    cmdProperAcc = cmdAcc + [0;0;g];

    % Saturate thrust magnitude (proper acceleration) and minimum vertical proper acc
    normProperAcc = norm(cmdProperAcc);
    if normProperAcc > params.maxProperAcc
        cmdProperAcc = cmdProperAcc * (params.maxProperAcc / normProperAcc);
        normProperAcc = params.maxProperAcc;
    end
    if cmdProperAcc(3) < params.minVerticalProperAcceleration
        cmdProperAcc(3) = params.minVerticalProperAcceleration;
        normProperAcc = norm(cmdProperAcc);
    end

    % Desired thrust direction
    thrustDir = cmdProperAcc / normProperAcc;

    % Build attitude with body z-axis = thrustDir, yaw = 0
    temp_up = [0; 0; 1];
    if abs(dot(thrustDir, temp_up)) > 0.99
        temp_up = [0; 1; 0];
    end
    x_axis = cross(temp_up, thrustDir); x_axis = x_axis / norm(x_axis);
    y_axis = cross(thrustDir, x_axis);
    R = [x_axis, y_axis, thrustDir];

    % Force in world frame (assumes body thrust along +z_b)
    T = m * normProperAcc;
    f = R * [0;0;T];
end

function p_des = desired_position(t)
    if t < 0.5
        p_des = [0; 0; 2];     % hover at (0,0)
    else
        p_des = [1; 2; 2];     % return to hover
    end
end
