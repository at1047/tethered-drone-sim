clear; clc; close all;

% --- Parameters ---
params.m_q = 0.7;
params.m_p = 0.7;
params.l = 0.5;
params.g = 9.81;

% --- Initial Conditions ---
% state s = [xq; yq; zq; dxq; dyq; dzq; xp; yp; zp; dxp; dyp; dzp]
xq0 = 0; yq0 = 0; zq0 = 2;
xp0 = xq0; yp0 = yq0; zp0 = zq0 - params.l; % Hanging straight down
s0 = [xq0; yq0; zq0; ...
      0; 0; 0; ...   % dxq, dyq, dzq
      xp0; yp0; zp0; ...
      0; 0; 0];      % dxp, dyp, dzp

tspan = [0 10];

% --- Run Simulation ---
[t, s] = ode45(@(t, s) drone_velocity_dynamics(t, s, params), tspan, s0);

% --- Plotting ---
figure('Name', 'Drone Velocity Control Simulation');
subplot(2,2,1);
plot(t, s(:, 1:3), 'LineWidth', 2);
legend('x_q', 'y_q', 'z_q');
ylabel('Drone Position (m)');
title('Drone Position'); grid on;

subplot(2,2,2);
plot(t, s(:, 7:9), 'LineWidth', 2);
legend('x_p', 'y_p', 'z_p');
ylabel('Payload Position (m)');
title('Payload Position'); grid on;

subplot(2,2,3);
% Calculate angles from redundant states for plotting
xq_plot = s(:, 1); yq_plot = s(:, 2); zq_plot = s(:, 3);
xp_plot = s(:, 7); yp_plot = s(:, 8); zp_plot = s(:, 9);
th_plot = real(acos( (zq_plot - zp_plot) ./ params.l ));
phi_plot = atan2(yp_plot - yq_plot, xp_plot - xq_plot);
plot(t, rad2deg([th_plot, phi_plot]), 'LineWidth', 2);
legend('\theta', '\phi');
ylabel('Angle (deg)');
title('Pendulum Angles (deg)'); grid on;

subplot(2,2,4);
plot3(s(:, 1), s(:, 2), s(:, 3), 'b-', 'LineWidth', 2); hold on;
plot3(s(:, 7), s(:, 8), s(:, 9), 'r--', 'LineWidth', 1.5);
legend('Drone', 'Payload');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Trajectory'); grid on; axis equal; view(3);

% --- Animation ---
dt_anim = 0.05;
t_anim = 0:dt_anim:t(end);
s_anim = interp1(t, s, t_anim);
animate(s_anim', dt_anim);

function ds = drone_velocity_dynamics(t, s, p)
    % Extract states (Convention: [pos; vel; pos; vel])
    xq = s(1); yq = s(2); zq = s(3);
    dxq = s(4); dyq = s(5); dzq = s(6);
    xp = s(7); yp = s(8); zp = s(9);
    dxp = s(10); dyp = s(11); dzp = s(12);

    % --- Command Input (Circular Trajectory) ---
    R = 1.0; omega = 1.5;
    vx_cmd = -R * omega * sin(omega * t);
    vy_cmd =  R * omega * cos(omega * t);
    vz_cmd = 0;
    
    ddxq_cmd = -R * omega^2 * cos(omega * t);
    ddyq_cmd = -R * omega^2 * sin(omega * t);
    ddzq_cmd = 0;

    % --- Physics (Using auto-generated redundant EOM) ---
    mq = p.m_q; mp = p.m_p; l = p.l; g = p.g;

    % --- Controller ---
    % Simple tracking controller for the drone to follow the velocity command
    kv = 20;
    f1 = mq * (ddxq_cmd + kv * (vx_cmd - dxq));
    f2 = mq * (ddyq_cmd + kv * (vy_cmd - dyq));
    % Hover force compensates for both masses in steady state
    f3 = mq * (ddzq_cmd + kv * (vz_cmd - dzq)) + (mq + mp) * g;

    % Evaluate redundant accelerations
    % (xq,yq,zq,xp,yp,zp,dxq,dyq,dzq,dxp,dyp,dzp,f1,f2,f3,m_q,m_p,l,g)
    ddq = auto_ddq_final(xq,yq,zq,xp,yp,zp,dxq,dyq,dzq,dxp,dyp,dzp,f1,f2,f3,mq,mp,l,g);

    % --- Assemble Derivative ---
    ds = [dxq;    % d/dt xq
          dyq;    % d/dt yq
          dzq;    % d/dt zq
          ddq(1); % d/dt dxq
          ddq(2); % d/dt dyq
          ddq(3); % d/dt dzq
          dxp;    % d/dt xp
          dyp;    % d/dt yp
          dzp;    % d/dt zp
          ddq(4); % d/dt dxp
          ddq(5); % d/dt dyp
          ddq(6)];% d/dt dzp

end
