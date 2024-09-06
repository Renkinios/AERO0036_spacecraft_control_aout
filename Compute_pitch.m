function [Iw_pitch,max_v,max_power] = Compute_pitch(Omega_max, beta, R, N, c, Iyy, degree, time, h, dx, graph, print_result)


color1 = '#FFD700';
color2 = '#FF6600';
color3 = '#04194E';
color4 = '#759AAB';
color5 = '#9E2A2B';


delta_theta = deg2rad(degree/2);
delta_t     = time/2;
t           = (0:h:time);

% Made in two phases:

p_dot = delta_theta * 2/delta_t^2;          % Scalar angular acceleration of the spacecraft (constant) [rad/s²], first phase
% fprintf("p_dot : %f \n",p_dot);
T     = Iyy * p_dot;                        % Scalar torque of the spacecraft (constant) [Nm], first phase
H_max = T * time/2;                         % Scalar maximum angular momentum of the Wheels/spacecraft [Nm·s], first phase


w_max_Rspeed = Omega_max * 2 * pi / 60;            % Maximum angular speed of the reaction wheel [rad/s]
Omega_dot    = Omega_max * 2 * pi / (60 * time/2); % Scalar wheel angular acceleration (constant) [rad/s²], first phase

% Inertia requirement for the wheel
Iw_pitch = T / (2 * Omega_dot * sin(beta)); % Scalar inertia of the wheel [kg·m²]

% Initial time
t_1 = (0:h:time/2);                % Time for the first phase [s]
t_2 = time/2+h:h:time;             % Time for the second phase [s] add time because >
t   = [t_1, t_2];                  % Total time for the roll [s]


% Satellite momentum profile [Nms]
H = [T*(t_1), 2*H_max - T*(t_2)] ;

if graph
    figure;
    plot(t, H, 'LineWidth', 2);
    xlabel('Time [s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des x
    ylabel('Angular momentum [Nm $\cdot$ s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des y
    grid on;
    
    hold off;
    
    % Définir le chemin du dossier et le créer s'il n'existe pas
    folder_path = 'figures/step_1_pitch';
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
    
    % Définir le chemin complet pour le fichier PDF
    save_path = fullfile(folder_path, 'angular_momentum_plot.pdf');
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Sauvegarder la figure en PDF
    print(gcf, save_path, '-dpdf', '-bestfit');
end


% Satellite torque profile [Nm]
T_vec = [ones(1, time/(2*h)) * T, 0, ones(1, time/(2*h)) * (-T)];

if graph
    figure;
    plot(t, T_vec, 'LineWidth', 2);
    xlabel('Time [s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des x
    ylabel('Torque [Nm]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des y
    grid on;
    hold off;
    
    % Définir le chemin du dossier et le créer s'il n'existe pas
    folder_path = 'figures/step_1_pitch';
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
    
   % Définir le chemin complet pour le fichier PDF
    save_path = fullfile(folder_path, 'torque_plot.pdf');
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Sauvegarder la figure en PDF
    print(gcf, save_path, '-dpdf', '-bestfit');
end


% Angular velocity
p1 = p_dot * t_1;                  % Satellite momentum profile during the first phase [Nm·s]
p2 = p1(end) + - p_dot * (t_2 - t_1(end)); % Satellite momentum profile during the second phase [Nm·s]
p  = [p1, p2];                     % Satellite momentum profile during the roll [Nm·s]

% Vector wheel angular velocity
% First wheel in this case is the 2nd and 4th which are going to be used.
omega2_t1 = Omega_dot * t_1;        % Wheel angular velocity during the first phase [rad/s]
omega2_t2 = omega2_t1(end) - Omega_dot * (t_2 - t_1(end)); % Wheel angular velocity during the second phase [rad/s]
Omega_2   = [omega2_t1, omega2_t2]; % Wheel angular velocity during the roll [rad/s]
Omega_dot_2v = [Omega_dot * ones(1,length(t_1)), -Omega_dot * ones(1,length(t_2))]; % Vector of the wheel angular acceleration [rad/s²]

% Second wheel used
omega4_t1 = -omega2_t1;             % Wheel angular velocity during the first phase [rad/s]
omega4_t2 = -omega2_t2;             % Wheel angular velocity during the second phase [rad/s]
Omega_4   = [omega4_t1, omega4_t2]; % Wheel angular velocity during the roll [rad/s]



% Electrical current profile : 
Qf2 = -c *( sin(beta) * p - Omega_2) ; 
i2  = (-Iw_pitch * Omega_dot_2v - Qf2) / N ;
e2  = R * i2 + N*(sin(beta) * p - Omega_2) ; 

% Compute max voltage

f_pitch      = @(N, R) (R*((- Iw_pitch * Omega_dot_2v - Qf2)/N) + N*(sin(beta)*p - Omega_2)); 
f_pitch_values = f_pitch(N, R);
[max_value, max_index] = max(f_pitch_values);
[min_value, min_index] = min(f_pitch_values);

index = 0 ;
if abs(max_value) < abs(min_value)
    index = min_index ;
else
    index = max_index ;
end
max_v = f_pitch_values(index);
voltage_circ = max_v - N*(sin(beta)*p(index) - Omega_2(index));
% fprintf(">>voltage_circ in roll : %f [kV] \n",abs(voltage_circ)/1000);

max_power = (voltage_circ)^2/R/10^6;




%======================%
        %Modif%
%======================%
i = T_vec / (2*N*sin(beta));
voltage = i * R;
voltage_t0 = voltage(1) / 1000;
%max current
i_max = max(abs(i));
voltage_max = i_max * R;
voltage_max = voltage_max / 1000; %  [kV]

power_max = voltage_max * i_max;
power_max = power_max / 1000; % [MW]


if graph 
    figure;
    hold on;
    grid on;
   
    plot(t, i, 'LineWidth', 2);
    ylabel('Electrical current [A]', 'Interpreter', 'latex');
    xlabel('Time [s]', 'Interpreter', 'latex');
    hold off;
    folder_path = 'figures/step_1_pitch';  % Path to the folder
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Create the folder if it does not exist
    end
    save_path = fullfile(folder_path, 'current_voltage_plot.pdf');  % Path to the file
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    print(gcf, save_path, '-dpdf', '-bestfit');
end


if print_result
    fprintf('>>-----------Pitch-----------<<\n')
    fprintf('>> T:                %f [Nm]    \n',T);
    fprintf('>> H_max:            %f [Nms]   \n',H_max);
    fprintf('>> w_max_Rspeed:     %f [rad/s] \n',w_max_Rspeed);
    fprintf('>> Iw:               %f [kgm^2] \n',Iw_pitch);
    fprintf('>> Current at t = 0: %f [A] \n', i(1));
    fprintf('>> Voltage at t = 0: %f [KV] \n', voltage_t0);
    fprintf('>> Max voltage: %f [KV] \n', voltage_max);
    fprintf('>> Max power: %f [MW] \n', power_max);
end
end
