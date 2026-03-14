close all; clear all; clc;

g = 9.81;

q0 = [0;
      0;
      2;];

dq0 = [0;
       0;
       0;];

qp0 = [0;
      0;
      1;];

dqp0 = [0;
       0;
       0;];

% Constants
params.mq = 0.7;
params.mp = 0.1;
params.L = 1;
params.e3 = [0;0;1];
params.g = -9.81;

% Controller & aero Constants
params.pos_natFreq = 2.0;                 % [rad/s] natural frequency for position loop
params.pos_damping = 0.7;                 % damping ratio
params.maxProperAcc = 13;                 % [m/s^2] max proper acceleration
params.minVerticalProperAcceleration = -1; % [m/s^2] min vertical proper acceleration
params.linDragCoeffB = [0.1286; 0.1286; 0.1286]; % body-frame linear drag [N/(m/s)]
params.payloadLinDragCoeffB = [0; 0; 0];    % payload linear drag [N/(m/s)] in world frame

% MPC parameters
params.mpc_horizon = 20;                  % Prediction horizon steps
params.mpc_dt = 0.1;                      % MPC time step [s]
params.mpc_use = true;                    % Toggle MPC on/off
params.compare_mode = true;               % Run both MPC on/off and compare
params.mpc_max_swing_angle = 15*pi/180;   % Max swing angle [rad] (hard constraint)
params.mpc_Q_pos = 50.0;                  % Cost weight: position tracking
params.mpc_Q_swing = 50.0;                % Cost weight: swing angle
params.mpc_Q_swing_rate = 50.0;           % Cost weight: swing angular rate
params.mpc_R = 1.0;                       % Cost weight: control effort
params.mpc_Qf_pos = 100.0;                % Terminal cost: position tracking
params.mpc_Qf_swing = 200.0;              % Terminal cost: swing angle (want zero)
params.mpc_Qf_swing_rate = 200.0;         % Terminal cost: swing rate (want zero)
params.mpc_Qf_vel = 100.0;                % Terminal cost: drone velocity (want zero)

% Initial conditions
rhat = [0;0;-1];

s_initial = [q0;dq0];
s0 = s_initial;

sp_initial = [qp0; dqp0];
sp0 = sp_initial;

dt = 0.002;
t_steps = 0:dt:10;

s_arr = zeros(6, length(t_steps));
s_arr(:, 1) = s0;

sp_arr = zeros(6, length(t_steps));
sp_arr(:, 1) = sp0;

T_arr = zeros(1, length(t_steps));

R_arr = zeros(3, length(t_steps));
R_arr(:,1) = rhat;

% MPC reference trajectory storage
q_ref_arr = zeros(3, length(t_steps));
dq_ref_arr = zeros(3, length(t_steps));
swing_angle_arr = zeros(1, length(t_steps));

% Storage for comparison mode (MPC off)
if params.compare_mode
    s_arr_no_mpc = zeros(6, length(t_steps));
    sp_arr_no_mpc = zeros(6, length(t_steps));
    swing_angle_arr_no_mpc = zeros(1, length(t_steps));
    s_arr_no_mpc(:, 1) = s0;
    sp_arr_no_mpc(:, 1) = sp0;
    
    % Storage for MPC without terminal cost
    s_arr_no_term = zeros(6, length(t_steps));
    sp_arr_no_term = zeros(6, length(t_steps));
    swing_angle_arr_no_term = zeros(1, length(t_steps));
    s_arr_no_term(:, 1) = s0;
    sp_arr_no_term(:, 1) = sp0;
end

% MPC runs at slower frequency than simulation
mpc_dt = params.mpc_dt;
mpc_update_counter = 0;
mpc_q_ref = q0;      % Current MPC reference position
mpc_dq_ref = dq0;    % Current MPC reference velocity

% Run comparison if enabled
if params.compare_mode
    fprintf('Running simulation with MPC ON...\n');
end

for i = 2:length(t_steps)
    t = t_steps(i);
    
    % MPC runs at lower frequency (every mpc_dt seconds)
    if params.mpc_use && (t >= mpc_update_counter * mpc_dt)
        % Get desired position from trajectory planner
        q_des = desired_position(t);
        
        % Run MPC to compute smooth reference trajectory
        [mpc_q_ref, mpc_dq_ref] = mpc_controller(params, t, s0, sp0, q_des);
        
        mpc_update_counter = mpc_update_counter + 1;
    end
    
    % PD controller tracks MPC reference (or desired if MPC off)
    if params.mpc_use
        q_target = mpc_q_ref;
        dq_target = mpc_dq_ref;
    else
        q_target = desired_position(t);
        dq_target = [0; 0; 0];
    end
    
    % Store MPC reference for plotting
    q_ref_arr(:, i) = q_target;
    dq_ref_arr(:, i) = dq_target;
    
    % Compute swing angle for monitoring
    rvec = sp0(1:3) - s0(1:3);
    swing_angle = atan2(norm(rvec(1:2)), -rvec(3));  % Angle from vertical
    swing_angle_arr(i) = swing_angle;
    
    % Simulate one step forward
    [s1, sp1, T1, rhat] = simulate(t, dt, s0, sp0, rhat, params, q_target, dq_target);
    s_arr(:, i) = s1;
    sp_arr(:, i) = sp1;
    T_arr(1, i) = T1;
    R_arr(:, i) = rhat;

    s0 = s1;
    sp0 = sp1;
end

% Run comparison simulation (MPC OFF)
if params.compare_mode
    fprintf('Running simulation with MPC OFF for comparison...\n');
    s0_no_mpc = s_initial;
    sp0_no_mpc = sp_initial;
    rhat_no_mpc = [0;0;-1];
    
    for i = 2:length(t_steps)
        t = t_steps(i);
        
        q_target_no_mpc = desired_position(t);
        dq_target_no_mpc = [0; 0; 0];
        
        [s1_no_mpc, sp1_no_mpc, ~, rhat_no_mpc] = simulate(t, dt, s0_no_mpc, sp0_no_mpc, rhat_no_mpc, params, q_target_no_mpc, dq_target_no_mpc);
        
        s_arr_no_mpc(:, i) = s1_no_mpc;
        sp_arr_no_mpc(:, i) = sp1_no_mpc;
        
        rvec_no_mpc = sp0_no_mpc(1:3) - s0_no_mpc(1:3);
        swing_angle_no_mpc = atan2(norm(rvec_no_mpc(1:2)), -rvec_no_mpc(3));
        swing_angle_arr_no_mpc(i) = swing_angle_no_mpc;
        
        s0_no_mpc = s1_no_mpc;
        sp0_no_mpc = sp1_no_mpc;
    end
    
    % Run MPC WITHOUT terminal cost for comparison
    fprintf('Running simulation with MPC ON (no terminal cost)...\n');
    params_no_term = params;
    params_no_term.mpc_Qf_pos = 0;
    params_no_term.mpc_Qf_swing = 0;
    params_no_term.mpc_Qf_swing_rate = 0;
    params_no_term.mpc_Qf_vel = 0;
    
    s0_no_term = s_initial;
    sp0_no_term = sp_initial;
    rhat_no_term = [0;0;-1];
    mpc_update_counter_no_term = 0;
    mpc_q_ref_no_term = q0;
    mpc_dq_ref_no_term = dq0;
    
    for i = 2:length(t_steps)
        t = t_steps(i);
        
        if t >= mpc_update_counter_no_term * mpc_dt
            q_des_nt = desired_position(t);
            [mpc_q_ref_no_term, mpc_dq_ref_no_term] = mpc_controller(params_no_term, t, s0_no_term, sp0_no_term, q_des_nt);
            mpc_update_counter_no_term = mpc_update_counter_no_term + 1;
        end
        
        [s1_no_term, sp1_no_term, ~, rhat_no_term] = simulate(t, dt, s0_no_term, sp0_no_term, rhat_no_term, params_no_term, mpc_q_ref_no_term, mpc_dq_ref_no_term);
        
        s_arr_no_term(:, i) = s1_no_term;
        sp_arr_no_term(:, i) = sp1_no_term;
        
        rvec_no_term = sp0_no_term(1:3) - s0_no_term(1:3);
        swing_angle_no_term = atan2(norm(rvec_no_term(1:2)), -rvec_no_term(3));
        swing_angle_arr_no_term(i) = swing_angle_no_term;
        
        s0_no_term = s1_no_term;
        sp0_no_term = sp1_no_term;
    end
    
    fprintf('Comparison complete!\n');
end

%% Plotting
figure(1); clf;
sgtitle('Drone Position')
subplot(3,1,1);
plot(t_steps, s_arr(1,:), 'r-', 'LineWidth', 2); hold on;
if params.mpc_use
    plot(t_steps, q_ref_arr(1,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'MPC ref');
end
xlabel('Time [s]'); ylabel('X [m]'); grid on; legend;

subplot(3,1,2);
plot(t_steps, s_arr(2,:), 'g-', 'LineWidth', 2); hold on;
if params.mpc_use
    plot(t_steps, q_ref_arr(2,:), 'g--', 'LineWidth', 1.5);
end
xlabel('Time [s]'); ylabel('Y [m]'); grid on;

subplot(3,1,3);
plot(t_steps, s_arr(3,:), 'b-', 'LineWidth', 2); hold on;
if params.mpc_use
    plot(t_steps, q_ref_arr(3,:), 'b--', 'LineWidth', 1.5);
end
xlabel('Time [s]'); ylabel('Z [m]'); grid on;

figure(2); clf;
sgtitle('Payload Position')
subplot(3,1,1); hold on;
if params.compare_mode
    plot(t_steps, sp_arr(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'MPC + Terminal');
    plot(t_steps, sp_arr_no_term(1,:), 'm-', 'LineWidth', 2, 'DisplayName', 'MPC (no terminal)');
    plot(t_steps, sp_arr_no_mpc(1,:), 'r--', 'LineWidth', 2, 'DisplayName', 'No MPC');
    legend('Location', 'best');
else
    plot(t_steps, sp_arr(1,:), 'r', 'LineWidth', 2);
end
xlabel('Time [s]'); ylabel('X [m]'); grid on;

subplot(3,1,2); hold on;
if params.compare_mode
    plot(t_steps, sp_arr(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'MPC + Terminal');
    plot(t_steps, sp_arr_no_term(2,:), 'm-', 'LineWidth', 2, 'DisplayName', 'MPC (no terminal)');
    plot(t_steps, sp_arr_no_mpc(2,:), 'r--', 'LineWidth', 2, 'DisplayName', 'No MPC');
    legend('Location', 'best');
else
    plot(t_steps, sp_arr(2,:), 'g', 'LineWidth', 2);
end
xlabel('Time [s]'); ylabel('Y [m]'); grid on;

subplot(3,1,3); hold on;
if params.compare_mode
    plot(t_steps, sp_arr(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'MPC + Terminal');
    plot(t_steps, sp_arr_no_term(3,:), 'm-', 'LineWidth', 2, 'DisplayName', 'MPC (no terminal)');
    plot(t_steps, sp_arr_no_mpc(3,:), 'r--', 'LineWidth', 2, 'DisplayName', 'No MPC');
    legend('Location', 'best');
else
    plot(t_steps, sp_arr(3,:), 'b', 'LineWidth', 2);
end
xlabel('Time [s]'); ylabel('Z [m]'); grid on;

figure(3); clf; hold on;
sgtitle('Payload Swing Angle - COMPARISON')
if params.compare_mode
    plot(t_steps, swing_angle_arr*180/pi, 'b-', 'LineWidth', 2, 'DisplayName', 'MPC + Terminal');
    plot(t_steps, swing_angle_arr_no_term*180/pi, 'm-', 'LineWidth', 2, 'DisplayName', 'MPC (no terminal)');
    plot(t_steps, swing_angle_arr_no_mpc*180/pi, 'r-', 'LineWidth', 2, 'DisplayName', 'No MPC');
    max_swing_mpc = max(swing_angle_arr)*180/pi;
    max_swing_no_term = max(swing_angle_arr_no_term)*180/pi;
    max_swing_no_mpc = max(swing_angle_arr_no_mpc)*180/pi;
    fprintf('\n=== SWING COMPARISON ===\n');
    fprintf('Max swing (MPC + terminal):    %.2f degrees\n', max_swing_mpc);
    fprintf('Max swing (MPC, no terminal):  %.2f degrees\n', max_swing_no_term);
    fprintf('Max swing (no MPC):            %.2f degrees\n', max_swing_no_mpc);
    fprintf('Terminal cost improvement over no terminal: %.1f%% reduction\n', (1 - max_swing_mpc/max_swing_no_term)*100);
    fprintf('MPC improvement over no MPC:               %.1f%% reduction\n', (1 - max_swing_mpc/max_swing_no_mpc)*100);
    fprintf('========================\n\n');
else
    plot(t_steps, swing_angle_arr*180/pi, 'b', 'LineWidth', 2, 'DisplayName', 'MPC + Terminal');
end
plot(t_steps, params.mpc_max_swing_angle*180/pi*ones(size(t_steps)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Max allowed');
xlabel('Time [s]'); ylabel('Swing Angle [deg]'); grid on; legend('Location', 'best');

figure(4); clf;
plot(t_steps(2:end), T_arr(1,2:end), 'r', 'LineWidth', 2)
title('Cable Tension'); xlabel('Time [s]'); ylabel('Tension [N]'); grid on;

figure(5); clf;
sgtitle('rhat over time')
subplot(3,1,1);
plot(t_steps, R_arr(1,:), 'g', 'LineWidth', 2)
xlabel('Time [s]'); ylabel('rhat_x'); grid on;
subplot(3,1,2);
plot(t_steps, R_arr(2,:), 'g', 'LineWidth', 2)
xlabel('Time [s]'); ylabel('rhat_y'); grid on;
subplot(3,1,3);
plot(t_steps, R_arr(3,:), 'g', 'LineWidth', 2)
xlabel('Time [s]'); ylabel('rhat_z'); grid on;

animate([s_arr;sp_arr], dt)

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [s1, sp1, T1, new_rhat] = simulate(t, dt, s0, sp0, rhat, params, q_ref, dq_ref)
% Simulate one timestep of the drone-payload system
% Uses full cable DAE dynamics with drag

    mq = params.mq;
    mp = params.mp;
    e3 = params.e3;
    g = params.g;
    L = params.L;
    
    q = s0(1:3);
    dq = s0(4:6);
    p = sp0(1:3);
    dp = sp0(4:6);
    
    % Keep cable direction valid
    rvec = p - q;
    if norm(rvec) > 1e-6
        rhat = rvec / norm(rvec);
    else
        rhat = [0; 0; -1];
    end
    
    [f, R_wb] = u_input(params, t, q, dq, q_ref, dq_ref);
    
    % Body-frame linear drag
    v_b = R_wb.' * dq;
    K = diag(params.linDragCoeffB);
    F_drag_b = -K * v_b;
    F_drag_w = R_wb * F_drag_b;
    
    % Payload drag in world frame
    Kp_payload = diag(params.payloadLinDragCoeffB);
    F_drag_p = -Kp_payload * dp;
    
    % DAE system: [mq*I, 0, -rhat; 0, mp*I, +rhat; -rhat', rhat', 0] * [ddq; ddp; T] = B
    A = [mq * eye(3), zeros([3,3]), -rhat;
         zeros([3,3]), mp * eye(3), +rhat;
         -rhat.', rhat.', 0];
    
    B = [ f + F_drag_w + e3 * mq * g;
          e3 * mp * g + F_drag_p;
          -norm(dp - dq)^2 / L ];
    
    % Robust solve
    if rcond(A) < 1e-10
        X = pinv(A) * B;
    else
        X = A \ B;
    end
    
    ddq = X(1:3);
    ddp = X(4:6);
    T1 = X(7);
    
    % Semi-implicit Euler integration
    new_dq = dq + dt*ddq;
    new_dp = dp + dt*ddp;
    new_q = q + dt*new_dq;
    new_p = p + dt*new_dp;
    
    rvec_new = new_p - new_q;
    if norm(rvec_new) > 1e-6
        new_rhat = rvec_new / norm(rvec_new);
    else
        new_rhat = [0; 0; -1];
    end
    
    sp1 = [new_p; new_dp];
    s1 = [new_q; new_dq];
end

function [f, R] = u_input(params, t, p, pdot, p_ref, pdot_ref)
% PD position controller
% Second-order (wn, zeta) position loop -> desired acceleration
% Saturates proper acceleration, constructs attitude

    m = params.mq + params.mp;
    g = 9.81;

    p_des = p_ref;
    pdot_des = pdot_ref;

    wn = params.pos_natFreq;
    zeta = params.pos_damping;

    Kp = wn^2 * eye(3);
    Kd = 2 * zeta * wn * eye(3);

    cmdAcc = Kp*(p_des - p) + Kd*(pdot_des - pdot);
    cmdProperAcc = cmdAcc + [0;0;g];

    % Saturate proper acceleration
    normProperAcc = norm(cmdProperAcc);
    if normProperAcc > params.maxProperAcc
        cmdProperAcc = cmdProperAcc * (params.maxProperAcc / normProperAcc);
        normProperAcc = params.maxProperAcc;
    end
    if cmdProperAcc(3) < params.minVerticalProperAcceleration
        cmdProperAcc(3) = params.minVerticalProperAcceleration;
        normProperAcc = norm(cmdProperAcc);
    end

    % Build attitude: body z-axis = thrust direction
    thrustDir = cmdProperAcc / normProperAcc;
    temp_up = [0; 0; 1];
    if abs(dot(thrustDir, temp_up)) > 0.99
        temp_up = [0; 1; 0];
    end
    x_axis = cross(temp_up, thrustDir); x_axis = x_axis / norm(x_axis);
    y_axis = cross(thrustDir, x_axis);
    R = [x_axis, y_axis, thrustDir];

    T = m * normProperAcc;
    f = R * [0;0;T];
end

function p_des = desired_position(t)
% Desired trajectory: step from (0,0,2) to (1,0,2) at t=0.5
    if t < 0.5
        p_des = [0; 0; 2];
    else
        p_des = [1; 0; 2];
    end
end

%% ========================================================================
%  MPC CONTROLLER - CONSTRAINED OPTIMIZATION VIA FMINCON
%  ========================================================================
function [q_ref, dq_ref] = mpc_controller(params, t, s0, sp0, q_des)
% Model Predictive Controller for payload swing suppression
%
% Solves the constrained optimization problem:
%
%   min_{u_0,...,u_{N-1}}  sum_{k=0}^{N-1} [ Q_pos*||q_k - q_des||^2
%                                           + Q_swing*theta_k^2
%                                           + Q_swing_rate*theta_dot_k^2
%                                           + R*||u_k||^2 ]
%                        + Qf_pos*||q_N - q_des||^2
%                        + Qf_swing*theta_N^2
%                        + Qf_swing_rate*theta_dot_N^2
%                        + Qf_vel*||dq_N||^2
%
%   subject to:
%       dynamics: (q_{k+1}, qp_{k+1}) = f(q_k, qp_k, u_k)  (via simulate)
%       theta_k <= theta_max   for all k = 0,...,N
%       |u_k| <= u_max         for all k = 0,...,N-1
%
% Decision variables: u_k = acceleration commands (3D) over N horizon steps
% The first control action is used to compute the reference for the PD controller.

    % Extract current states
    q = s0(1:3);
    dq = s0(4:6);
    qp = sp0(1:3);
    dqp = sp0(4:6);
    
    % MPC parameters
    N = params.mpc_horizon;
    dt_mpc = params.mpc_dt;
    max_swing = params.mpc_max_swing_angle;
    
    % Stage cost weights
    Q_pos = params.mpc_Q_pos;
    Q_swing = params.mpc_Q_swing;
    Q_swing_rate = params.mpc_Q_swing_rate;
    R_weight = params.mpc_R;
    
    % Terminal cost weights
    Qf_pos = params.mpc_Qf_pos;
    Qf_swing = params.mpc_Qf_swing;
    Qf_swing_rate = params.mpc_Qf_swing_rate;
    Qf_vel = params.mpc_Qf_vel;
    
    % Decision variables: u = [u_0; u_1; ...; u_{N-1}], each u_k in R^3
    % Total: 3*N variables
    u_init = zeros(3*N, 1);
    
    % Bounds on control (acceleration limits)
    max_accel = params.maxProperAcc;
    u_lb = -max_accel * ones(3*N, 1);
    u_ub =  max_accel * ones(3*N, 1);
    
    % --------------------------------------------------------------------
    % Nested cost function (accesses parent scope variables)
    % --------------------------------------------------------------------
    function J = cost_fn(u_flat)
        u_seq = reshape(u_flat, 3, N);
        
        % Forward simulate the full dynamics over the horizon
        [q_traj, dq_pred, qp_traj, ~, theta_traj] = predict_trajectory(u_seq, ...
            q, dq, qp, dqp, t, N, dt_mpc, params);
        
        J = 0;
        
        % Stage costs (k = 1 to N, i.e. steps 0 to N-1)
        for k = 1:N
            % Position tracking
            e_pos = q_traj(:, k) - q_des;
            J = J + Q_pos * (e_pos' * e_pos);
            
            % Swing angle penalty
            J = J + Q_swing * theta_traj(k)^2;
            
            % Swing rate penalty
            if k > 1
                dtheta = (theta_traj(k) - theta_traj(k-1)) / dt_mpc;
                J = J + Q_swing_rate * dtheta^2;
            end
            
            % Control effort cost
            J = J + R_weight * (u_seq(:,k)' * u_seq(:,k));
        end
        
        % Terminal cost at step N (want: at goal, no swing, no motion)
        e_pos_f = q_traj(:, N+1) - q_des;
        J = J + Qf_pos * (e_pos_f' * e_pos_f);
        J = J + Qf_swing * theta_traj(N+1)^2;
        dtheta_f = (theta_traj(N+1) - theta_traj(N)) / dt_mpc;
        J = J + Qf_swing_rate * dtheta_f^2;
        J = J + Qf_vel * (dq_pred(:, N+1)' * dq_pred(:, N+1));
    end
    
    % --------------------------------------------------------------------
    % Nested constraint function (accesses parent scope variables)
    % --------------------------------------------------------------------
    function [c, ceq] = constraint_fn(u_flat)
        u_seq = reshape(u_flat, 3, N);
        
        % Forward simulate
        [q_traj, ~, qp_traj, ~, theta_traj] = predict_trajectory(u_seq, ...
            q, dq, qp, dqp, t, N, dt_mpc, params);
        
        % Inequality constraint: theta_k <= max_swing  =>  theta_k - max_swing <= 0
        c = theta_traj(:) - max_swing;
        
        % No equality constraints
        ceq = [];
    end
    
    % --------------------------------------------------------------------
    % Solve with fmincon (SQP)
    % --------------------------------------------------------------------
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-4, ...
        'StepTolerance', 1e-6);
    
    try
        u_opt = fmincon(@cost_fn, u_init, [], [], [], [], u_lb, u_ub, ...
                        @constraint_fn, options);
        
        % Extract first control action and compute reference for PD
        u0 = u_opt(1:3);
        q_ref = q + dt_mpc * dq + 0.5 * dt_mpc^2 * u0;
        dq_ref = dq + dt_mpc * u0;
        
    catch ME
        % Fallback: just track desired position directly
        warning('MPC optimization failed: %s. Using fallback.', ME.message);
        q_ref = q_des;
        dq_ref = [0; 0; 0];
    end
    
    % Safety check
    if any(isnan(q_ref)) || any(isnan(dq_ref))
        q_ref = q;
        dq_ref = dq;
    end
end

%% ========================================================================
%  TRAJECTORY PREDICTION (used by MPC cost and constraints)
%  ========================================================================
function [q_traj, dq_traj, qp_traj, dqp_traj, theta_traj] = ...
    predict_trajectory(u_seq, q0, dq0, qp0, dqp0, t0, N, dt_mpc, params)
% Forward-simulate the drone-payload system over N steps using the full
% dynamics (calls simulate at each step).
%
% Inputs:
%   u_seq    - 3 x N control sequence (acceleration commands)
%   q0, dq0  - initial drone position/velocity
%   qp0, dqp0 - initial payload position/velocity
%   t0       - current time
%   N        - horizon length
%   dt_mpc   - MPC time step
%   params   - system parameters
%
% Outputs:
%   q_traj   - 3 x (N+1) predicted drone positions
%   dq_traj  - 3 x (N+1) predicted drone velocities
%   qp_traj  - 3 x (N+1) predicted payload positions
%   dqp_traj - 3 x (N+1) predicted payload velocities
%   theta_traj - 1 x (N+1) predicted swing angles

    q_traj = zeros(3, N+1);
    dq_traj = zeros(3, N+1);
    qp_traj = zeros(3, N+1);
    dqp_traj = zeros(3, N+1);
    theta_traj = zeros(1, N+1);
    
    q_traj(:,1) = q0;
    dq_traj(:,1) = dq0;
    qp_traj(:,1) = qp0;
    dqp_traj(:,1) = dqp0;
    
    % Initial swing angle
    rvec = qp0 - q0;
    if norm(rvec) > 1e-6
        theta_traj(1) = atan2(norm(rvec(1:2)), -rvec(3));
        rhat_local = rvec / norm(rvec);
    else
        theta_traj(1) = 0;
        rhat_local = [0; 0; -1];
    end
    
    for k = 1:N
        q_k = q_traj(:, k);
        dq_k = dq_traj(:, k);
        u_k = u_seq(:, k);
        
        % The MPC control u_k is an acceleration command.
        % Convert to a reference position/velocity that the PD controller tracks.
        q_ref_k = q_k + dt_mpc * dq_k + 0.5 * dt_mpc^2 * u_k;
        dq_ref_k = dq_k + dt_mpc * u_k;
        
        % Simulate one MPC step using full dynamics
        s0_k = [q_traj(:,k); dq_traj(:,k)];
        sp0_k = [qp_traj(:,k); dqp_traj(:,k)];
        t_k = t0 + (k-1) * dt_mpc;
        
        [s1_k, sp1_k, ~, rhat_local] = simulate(t_k, dt_mpc, s0_k, sp0_k, ...
                                                  rhat_local, params, q_ref_k, dq_ref_k);
        
        q_traj(:, k+1) = s1_k(1:3);
        dq_traj(:, k+1) = s1_k(4:6);
        qp_traj(:, k+1) = sp1_k(1:3);
        dqp_traj(:, k+1) = sp1_k(4:6);
        
        % Compute swing angle
        rvec_k = qp_traj(:, k+1) - q_traj(:, k+1);
        if norm(rvec_k) > 1e-6
            theta_traj(k+1) = atan2(norm(rvec_k(1:2)), -rvec_k(3));
        else
            theta_traj(k+1) = 0;
        end
    end
end