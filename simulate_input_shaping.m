clear
clc
close all

q0 = [0;
      deg2rad(15);] ;
dq0 = [0;
       0] ;
x0 = [q0;
      dq0] ;

m_c = 0.10;
m_l = 0.01;
l = 1;
g = 9.81;

wn = sqrt(g / l); 
Tn = 2 * pi / wn;      % Natural Period

shaper_times = [0, Tn/2];
shaper_amps  = [0.5, 0.5];

fprintf('Shaper Impulse Times: %.3f s and %.3f s\n', shaper_times(1), shaper_times(2));

dt = 0.005;
t_final = 10;
time_vec = 0:dt:t_final;

u_raw = zeros(size(time_vec));
u_raw(time_vec >= 0.5 & time_vec < 2.5) = 2.0; 

u_shaped = zeros(size(time_vec));

for i = 1:length(shaper_times)
    idx_shift = round(shaper_times(i) / dt);
    
    u_shifted = [zeros(1, idx_shift), u_raw(1:end-idx_shift)];
    u_shaped = u_shaped + (u_shifted * shaper_amps(i));
end

% Simulate Unshaped
[t1, y1] = ode45(@(t,y) dynamics(t, y, time_vec, u_raw), [0 t_final], x0);

% Simulate Shaped
[t2, y2] = ode45(@(t,y) dynamics(t, y, time_vec, u_shaped), [0 t_final], x0);

figure('Color','w', 'Position', [100 100 800 600]);

% Plot Inputs
subplot(3,1,1);
plot(time_vec, u_raw, 'r--', 'DisplayName', 'Unshaped', 'LineWidth', 1.5); hold on;
plot(time_vec, u_shaped, 'b-', 'DisplayName', 'ZV Shaped', 'LineWidth', 2);
ylabel('Force (N)');
title('Input Shaping: Breaking the Step into Two');
grid on;

% Plot Cart Position
subplot(3,1,2);
plot(t1, y1(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t2, y2(:,1), 'b-', 'LineWidth', 2);
ylabel('Cart Position (m)');
grid on;

% Plot Pendulum Swing (The critical part!)
subplot(3,1,3);
plot(t1, rad2deg(y1(:,2)), 'r--', 'DisplayName', 'Unshaped', 'LineWidth', 1.5); hold on;
plot(t2, rad2deg(y2(:,2)), 'b-', 'DisplayName', 'ZV Shaped', 'LineWidth', 2);
ylabel('Swing Angle (deg)'); xlabel('Time (s)');
grid on;

% Animate the system using the computed solution
% animate_system(t2, y2, 1);

function dx = dynamics(t, s, t_vec, u_vec)
    q = s(1:2) ; dq = s(3:4) ;
    x = q(1) ; th = q(2) ;
    dx = dq(1) ; dth = dq(2) ;

    u = interp1(t_vec, u_vec, t, 'linear', 0);

    % System Constants
    m_c = 0.10;
    m_l = 0.01;
    l = 1;
    g = 9.81;

    D = auto_D(l,m_c,m_l,th);
    C = auto_C(dth,l,m_l,th);
    G = auto_G(g,l,m_l,th);
    % B = auto_B() ;
    B = [1; 0];

    f = [dq ;
         inv(D)*(-C*dq - G)] ;
    g = [zeros(2,1); inv(D)*B] ;

    dx = f + g*u ;
end
