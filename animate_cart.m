
function animate_cart(t, state, l)
    % Inputs:
    % t: time vector
    % state: state matrix [x, theta, dx, dtheta]
    % l: length of pendulum rod
    
    % 1. Setup Figure
    figure('Color', 'w', 'Name', 'Cart-Pendulum Animation');
    
    % Visual Parameters
    W = 0.4;  % Cart Width
    H = 0.2;  % Cart Height
    mr = 0.05; % Mass Radius (for the bob)
    
    % Create Plot Objects (Initialize at t=0)
    % Ground line
    plot([-10 10], [-H/2 -H/2], 'k-', 'LineWidth', 2); hold on;
    
    % Cart (Rectangle)
    % Defined by 4 corners: [x_left, x_right, x_right, x_left]
    h_cart = patch('XData', [0 0 0 0], 'YData', [0 0 0 0], ...
                   'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'k');
    
    % Rod (Line)
    h_rod = plot([0 0], [0 0], 'k-', 'LineWidth', 2);
    
    % Bob (Circle/Marker)
    h_bob = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Time Text
    h_text = text(0, 0, '', 'FontSize', 12);

    % Axis Settings
    axis equal; grid on;
    % Fix axis limits so the camera follows the cart but doesn't zoom in/out
    ylim([-l*1.5, l*1.5]); 
    xlabel('Position (m)'); ylabel('Height (m)');
    
    % 2. Animation Loop
    % We skip frames if data is too dense to keep animation fast
    stride = max(1, floor(length(t) / 300)); 
    
    for i = 1:stride:length(t)
        % --- Dynamics Kinematics ---
        x_c = state(i, 1);      % Cart Position
        th  = state(i, 2);      % Angle
        
        % Pendulum Bob Position (Relative to Cart)
        x_p = x_c + l*sin(th);
        y_p = -l*cos(th);       % Down is negative Y
        
        % --- Update Graphics ---
        
        % Update Cart (Centered at x_c)
        set(h_cart, 'XData', [x_c-W/2, x_c+W/2, x_c+W/2, x_c-W/2], ...
                    'YData', [-H/2, -H/2, H/2, H/2]);
        
        % Update Rod (From Cart center to Bob)
        set(h_rod, 'XData', [x_c, x_p], 'YData', [0, y_p]);
        
        % Update Bob
        set(h_bob, 'XData', x_p, 'YData', y_p);
        
        % Update Camera/Text
        set(h_text, 'String', sprintf('Time: %.2fs', t(i)), ...
                    'Position', [x_c - 1, l + 0.2]);
        
        % Follow the cart with the camera
        xlim([x_c - 2, x_c + 2]);
        
        % Force Draw
        drawnow; 
        
        % Optional: Pause to match real-time
        % pause(0.01); 
    end
end