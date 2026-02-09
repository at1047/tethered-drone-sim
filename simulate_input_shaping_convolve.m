clear; clc; close all;

% Declare Globals

% Initial Conditions
q0 = [0; deg2rad(0)];
dq0 = [0; 0];
x0 = [q0; dq0];

% Parameters
params.m_c = 0.7;
params.m_l = 0.1;
params.l = 0.5;
params.g = 9.81;
params.u_des = 3;

% Shaper Calculations
params.wn = sqrt(params.g / params.l); 
params.Tn = 2 * pi / params.wn; 
dT = params.Tn / 2;

% Simulation Setup
params.dt = 1/500;
t_final = 10;
t_vec = 0:params.dt:t_final;

% Initialize Global Arrays
u_arr = zeros(size(t_vec));
u_convolve = zeros(size(t_vec));

% Simulate
t = 0;
x = x0;
y1 = [];
t1 = [];

periodOffBoardMainLoop = 1/50;
timerOffBoardMainLoop = 0;

while t < t_final
    
    idx = round(t / params.dt) + 1;

    % if timerOffBoardMainLoop > periodOffBoardMainLoop 
    %     timerOffBoardMainLoop = 0;
    % 
    %     idx_shift = round((params.Tn/2) / params.dt);
    % 
    %     u_raw = u_input(t, params, x(1:2), x(3:4));
    % 
    %     % u_arr(idx:idx_shift) = u_raw;
    %     % Impulse 1 (Half energy now)
    %     idxEnd_1 = min(idx + periodOffBoardMainLoop/params.dt, length(u_convolve));
    %     u_convolve(idx:idxEnd_1) = u_convolve(idx:idxEnd_1) + u_raw/2;
    % 
    %     % Impulse 2 (Half energy delayed by Tn/2)
    % 
    %     if (idx + idx_shift) <= length(u_convolve)
    %         idxEnd_2 = min(idx + idx_shift + periodOffBoardMainLoop/params.dt, length(u_convolve));
    %         u_convolve(idx + idx_shift:idxEnd_2) = u_convolve(idx + idx_shift:idxEnd_2) + u_raw/2;
    %     end
    % end

    if timerOffBoardMainLoop > periodOffBoardMainLoop 
        timerOffBoardMainLoop = 0;

        idx_shift = round((params.Tn/2) / params.dt);

        u_raw = u_input(t, params, x(1:2), x(3:4));

        % u_arr(idx:idx_shift) = u_raw;
        % Impulse 1 (Half energy now)
        idxEnd_1 = min(idx + periodOffBoardMainLoop/params.dt, length(u_convolve));
        u_convolve(idx:idxEnd_1) = u_convolve(idx:idxEnd_1) + u_raw/2;
    
        % Impulse 2 (Half energy delayed by Tn/2)
        
        if (idx + idx_shift) <= length(u_convolve)
            idxEnd_2 = min(idx + idx_shift + periodOffBoardMainLoop/params.dt, length(u_convolve));
            u_convolve(idx + idx_shift:idxEnd_2) = u_convolve(idx + idx_shift:idxEnd_2) + u_raw/2;
        end
    end

    dx = dynamics(t, x, u_convolve(idx), params);
    x = x + dx*params.dt;
    t1 = [t1; t];
    y1 = [y1; x.'];

    t = t + params.dt;
    timerOffBoardMainLoop = timerOffBoardMainLoop + params.dt;
end

% --- Plotting ---
figure('Color','w', 'Position', [100 100 800 800]);

subplot(3,1,1);
plot(t1, y1(:,1), 'b-', 'LineWidth', 2);
ylabel('Cart Position (m)'); grid on;

subplot(3,1,2);
plot(t1, rad2deg(y1(:,2)), 'r-', 'LineWidth', 2);
ylabel('Swing Angle (deg)'); grid on;

subplot(3,1,3);
% Plotting the global arrays against the time vector
plot(t_vec, u_arr, 'k--', 'DisplayName', 'Raw Command'); hold on;
plot(t_vec, u_convolve, 'b-', 'LineWidth', 1.5, 'DisplayName', 'ZV Shaped');
ylabel('Control Force (N)'); xlabel('Time (s)');
legend; grid on;

function dx = dynamics(t, s, u, params)
    
    q = s(1:2); dq = s(3:4);
    th = q(2); dth = q(2);

    % Physics
    D = auto_D(params.l, params.m_c, params.m_l, th);
    C = auto_C(dth, params.l, params.m_l, th);
    G = auto_G(params.g, params.l, params.m_l, th);
    B = [1; 0];

    f = [dq; inv(D)*(-C*dq - G)];
    g = [zeros(2,1); inv(D)*B];
    dx = f + g*u;
end

function u = u_input(t, params, q, dq)
    x = q(1); dx = dq(1);
    x_des = 3; 
    kp = 1; kd = 0.25;
    u = -kp * (x - x_des) - kd * (dx - 0);
end