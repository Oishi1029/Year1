% Define phasor data 
% Define phasor data 
voltage_mag = 120;       % Magnitude of the voltage phasor
voltage_angle = 0;       % Angle of the voltage phasor (degrees)
voltage_uncertainty = 5; % Uncertainty in voltage magnitude
current_mag = 5;         % Magnitude of the current phasor
current_angle = -30;     % Angle of the current phasor (lagging)
current_uncertainty = 0.2;% Uncertainty in current magnitude

% Convert angles to radians
voltage_angle_rad = deg2rad(voltage_angle);
current_angle_rad = deg2rad(current_angle);

% Calculate phasor endpoints (central values)
vx = voltage_mag * cos(voltage_angle_rad);
vy = voltage_mag * sin(voltage_angle_rad);
ix = current_mag * cos(current_angle_rad);
iy = current_mag * sin(current_angle_rad);

% Create the figure
figure;
hold on;  

% Voltage uncertainty phasors
for offset = [-voltage_uncertainty, voltage_uncertainty]
    new_vx = vx + offset;
    new_vy = vy; % No change in y for voltage uncertainty 
    quiver(0, 0, new_vx, new_vy, 'LineWidth', 1, 'Color', 'b', 'LineStyle', '--');
end

% Current uncertainty phasors
for offset = [-current_uncertainty, current_uncertainty]
    new_ix = ix + offset;
    new_iy = iy; % No change in y for current uncertainty 
    quiver(0, 0, new_ix, new_iy, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
end

% ... (Your existing code to define phasors and plot them) ...

% Labels for voltage phasors
text(vx, vy, 'V', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Central voltage  
for offset = [-voltage_uncertainty, voltage_uncertainty]
    new_vx = vx + offset;
    text(new_vx,  new_vy, 'V (?)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Uncertain phasors
end

% Labels for current phasors
text(ix, iy, 'I', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Central current
for offset = [-current_uncertainty, current_uncertainty]
    new_ix = ix + offset;
    text(new_ix, new_iy, 'I (?)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Uncertain phasors
end

% ... (Labels, title, legend, etc.) ...

% Plot the main phasors with thicker lines
quiver(0, 0, vx, vy, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Voltage');
quiver(0, 0, ix, iy, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Current');

% ... (Labels, title, legend, etc.) ...

% ... (Rest of your code for labels, title, legend, etc.) ...

% Labels, title, and legend 
xlabel('Real');
ylabel('Imaginary');
title('Phasor Diagram with Uncertainties');
legend();  
grid on; 
axis equal; % Ensures correct aspect ratio for phasors
