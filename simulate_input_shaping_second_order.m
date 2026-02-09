clear
clc
close all

q0 = [0;
      deg2rad(0);] ;
dq0 = [0;
       0] ;
x0 = [q0;
      dq0] ;

params.m_c = 0.7;
params.m_l = 0.1;
params.l = 0.5;
params.g = 9.81;

params.wn = sqrt(params.g / params.l); 
params.Tn = 2 * pi / params.wn;      % Natural Period

params.u_des = 3;

shaper_times = [0, params.Tn/2];
shaper_amps  = [0.5, 0.5];

fprintf('Shaper Impulse Times: %.3f s and %.3f s\n', shaper_times(1), shaper_times(2));

params.dt = 0.005;
t_final = 20;
time_vec = 0:params.dt:t_final;

% params.shaper_times = shaper_times;


% Simulate Unshaped
% [t1, y1] = ode45(@(t,y) dynamics(t, y, params), [0 t_final], x0);

% Simulate Shaped
params.t_arr = time_vec;
params.u_convolve = zeros(size(time_vec));
params.u_arr = zeros(size(time_vec));
[t2, y2] = ode45(@(t,y) dynamics(t, y, params), [0 t_final], x0);

figure('Color','w', 'Position', [100 100 800 600]);

% Plot Inputs
% subplot(3,1,1);
% plot(time_vec, u_raw, 'r--', 'DisplayName', 'Unshaped', 'LineWidth', 1.5); hold on;
% plot(time_vec, u_shaped, 'b-', 'DisplayName', 'ZV Shaped', 'LineWidth', 2);
% ylabel('Force (N)');
% title('Input Shaping: Breaking the Step into Two');
% grid on;

% Plot Cart Position
subplot(3,1,1);
% plot(t1, y1(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t2, y2(:,1), 'b-', 'LineWidth', 2);
ylabel('Cart Position (m)');
grid on;

% Plot Pendulum Swing (The critical part!)
subplot(3,1,2);
% plot(t1, rad2deg(y1(:,2)), 'r--', 'DisplayName', 'Unshaped', 'LineWidth', 1.5); hold on;
plot(t2, rad2deg(y2(:,2)), 'b-', 'DisplayName', 'ZV Shaped', 'LineWidth', 2);
ylabel('Swing Angle (deg)'); xlabel('Time (s)');
grid on;

subplot(3,1,3);
% plot(t1, rad2deg(y1(:,2)), 'r--', 'DisplayName', 'Unshaped', 'LineWidth', 1.5); hold on;
plot(params.t_arr, params.u_arr, 'r--', 'DisplayName', 'raw', 'LineWidth', 2);
plot(params.t_arr, params.u_convolve, 'b-', 'DisplayName', 'ZV Shaped', 'LineWidth', 2);
ylabel('Swing Angle (deg)'); xlabel('Time (s)');
grid on;

% Animate the system using the computed solution
% animate_system(t2, y2, 1);

function dx = dynamics(t, s, params)
    q = s(1:2) ; dq = s(3:4) ;
    x = q(1) ; th = q(2) ;
    dx = dq(1) ; dth = dq(2) ;

    % u = interp1(t_vec, u_vec, t, 'linear', 0);

    D = auto_D(params.l,params.m_c,params.m_l,th);
    C = auto_C(dth,params.l,params.m_l,th);
    G = auto_G(params.g,params.l,params.m_l,th);
    % B = auto_B() ;
    B = [1; 0];

    idx = round(t / params.dt)+1;
    idx_shift = round((params.Tn/2) / params.dt) + 1;

    u_raw = u_input(t, params, q, dq);

    params.u_arr(idx) = u_raw;
    params.u_convolve(idx) = params.u_convolve(idx) + u_raw/2;
    params.u_convolve(idx_shift) = params.u_convolve(idx_shift) + u_raw/2;

    u = interp1(params.t_arr, params.u_convolve, t, 'linear', 0);

    f = [dq ;
         inv(D)*(-C*dq - G)] ;
    g = [zeros(2,1); inv(D)*B] ;

    dx = f + g*u ;
end

function u = u_input(t, params, q, dq)

    x = q(1); dx = dq(1);

    % tau = params.Tn/2;

    % if t < tau
    % x_des = 1.5; dx_des = 0;
    % else
    % x_des = 3; dx_des = 0;
    % end

    x_des = 3; dx_des = 0;

    kp = 1;
    kd = 0.25;
    
    u = -kp * (x-x_des) -kd * (dx-dx_des);

end