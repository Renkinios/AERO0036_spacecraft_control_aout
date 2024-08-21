function [Iw_roll,max_v,max_p] = Compute_roll(Ixx,beta,R,N,c,time_step,Omega_max,do_fig,print_result)

t_roll     = (0:time_step:10) ;      % Time to be able to change its orientation of 90 degrees in roll (+x) [s]
t_r        = 10;                     % Total time for rotation in roll [s]
dx_roll    = 5 ;                     % First interval time derivate
delta_roll = deg2rad(45);            % Rotation angle in roll [deg]
delta_t    = dx_roll;                % Time interval to perform the rotation [s] 


p_dot        = delta_roll*2/delta_t^2;         % Scalar .Angular acceleration of the spacecraft (cst) [rad s-2], first phase by kinetic equation
% fprintf("p_dot : %f \n",p_dot);
T            = Ixx*p_dot;                      % Scalar .Torque of the spacecraft (cst) [Nm], first phase by euler's equation
H_max        = T*dx_roll;                      % Scalar .Maximum angular momentum of the Wheels/spacecraft [Nm s], first phase

fprintf("H_max : %f \n",H_max);
w_max_Rspeed = 7000*2*pi/60;                   % Scalar .Maximum rotation speed of the wheel [rad s-1], first phase

% initial time
t_1 = 0:time_step:5;              % Time for the first phase [s]
t_2 = 5 + time_step:time_step:10; % Time for the second phase [s] add time becausse >
t   = [t_1, t_2];                 % Total time for the roll [s]

H_start = T * t_1;                     % Satellite momentum progile during the first phase [Nm s]
H_end   = 2 * H_max - T *(t_2);        % Satellite momentum progile during the second phase [Nm s]
H       = [H_start, H_end];            % Satellite momentum progile during the roll [Nm s]

if do_fig 
    figure;
    plot(t, H, 'LineWidth', 2);
    xlabel('Time [s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des x
    ylabel('Angular momentum [Nm $\cdot$ s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des y
    grid on;
    
    % Création du dossier si nécessaire
    folder_path = 'figures/step_1_roll';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end
    
    % Sauvegarde de la figure en PDF
    save_path = fullfile(folder_path, 'angular_momentum_plot.pdf');  % Chemin complet du fichier de sauvegarde
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Sauvegarder la figure en PDF
    print(gcf, save_path, '-dpdf', '-bestfit');
end

% Torque of the spacecraft
T_vec = [ones(1, length(t_1)) * T, ones(1, length(t_2)) * (-T)];  % Vector .Torque of the spacecraft [Nm]

if do_fig 
    figure;
    plot(t, T_vec, 'LineWidth', 2);
    xlabel('Time [s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des x
    ylabel('Torque [Nm]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des y
    grid on;
    
    % Création du dossier si nécessaire
    folder_path = 'figures/step_1_roll';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
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

Omega_dot    = Omega_max * 2 * pi/60 / dx_roll;    % Scalar .Wheel angular acceleration (cst) [rad s-2], first phase
Iw_roll      = T/(2*Omega_dot*sin(beta));          % Scalar .Inertia of the wheel [kg m²], first phase
% ? 

if print_result
    fprintf('>>------------Roll-----------<<\n');
    fprintf('>> T:                %f [Nm] \n',T);
    fprintf('>> H_max:            %f [Nms] \n',H_max);
    fprintf('>> w_max_Rspeed:     %f [rad/s] \n',w_max_Rspeed);
    fprintf('>> Iw:               %f [kgm^2] \n',Iw_roll);
end

% Angulaire acceleration of the wheels
alpha1    = Omega_dot;                                                   % Scalar .Wheel angular acceleration [rad s-2], first phase
alpha2    = - Omega_dot;                                                 % Scalar .Wheel angular acceleration [rad s-2], first phase
alpha     = [alpha1*ones(1,length(t_1)), alpha2*ones(1,length(t_2))]; 

% Angular velocity 
p_1 = p_dot * t_1;                                    % Vector .Angular velocity of the spacecraft during the first phase [rad s-1]
p_2 = p_1(end) - p_dot * (t_2 - t_1(end));            % Vector .Angular velocity of the spacecraft during the second phase [rad s-1]
p   = [p_1, p_2];                                     % Vector .Angular velocity of the spacecraft during the roll [rad s-1]

% first weel 
Omega1_t1  = Omega_dot * t_1;                                % Vector .Angular velocity of the wheel 1 during the first phase [rad s-1]
Omega1_t2  = Omega1_t1(end) - Omega_dot * (t_2 - t_1(end));  % Vector .Angular velocity of the wheel 1 during the second phase [rad s-1]
Omega1     = [Omega1_t1, Omega1_t2];  


if do_fig
    figure ;
    plot(t, Omega1, 'LineWidth', 2)
    xlabel('Time [s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des x
    ylabel('Angular velocity [rad/s]', 'Interpreter', 'latex');  % Utilisation de LaTeX pour l'axe des y
    grid on;
    % Création du dossier si nécessaire
    folder_path = 'figures/step_1_roll';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end

    % Sauvegarde de la figure en PDF
    save_path = fullfile(folder_path, 'angular_velocity_wheel_1_plot.pdf');  % Chemin complet du fichier de sauvegarde
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Sauvegarder la figure en PDF
    print(gcf, save_path, '-dpdf', '-bestfit');
end
% second weel

Omega2_t1 = - Omega1_t1;                              % Vector .Angular velocity of the wheel 2 during the first phase [rad s-1]
Omega2_t2 = - Omega1_t2;                              % Vector .Angular velocity of the wheel 2 during the second phase [rad s-1]
Omega2    = [Omega2_t1, Omega2_t2];                   % Vector .Angular velocity of the wheel 2 during the roll [rad s-1]
% Electrical current profile


Qf1 =  - c*(sin(beta)*p - Omega1);                    % Vector .Friction torque [Nm] using Nexton at spacecraft

% Créer le plot
if do_fig
    figure();
    hold on;
    grid on;

    % Premier axe y pour la tension
    yyaxis left;
    plot(t,Iw_roll * alpha, 'LineWidth', 2);
    ylabel('p [V]', 'Interpreter', 'latex');

    % Deuxième axe y pour la force contre-EMF
    yyaxis right;
    plot(t, Qf1, 'LineWidth', 2);
    ylabel('Qf1 [A]', 'Interpreter', 'latex');

    % Axe x pour le temps
    xlabel('Time [s]', 'Interpreter', 'latex');

    % Légendes
    legend('Tension (p)', 'Force contre-EMF (Qf1)', 'Location', 'best');

    hold off;
end

% Obtenir la taille
[rows, cols] = size(Qf1);

i1  = (-Iw_roll * alpha - Qf1)/N;                     % Vector .Electrical current [A] newton as reaction wheel
e1  = R*i1 + N*(sin(beta) * p - Omega1);              % Vector .Armature volrage (control input) [V]


% Computation max voltage 
f_roll = @(N, R) (R * ((- Iw_roll * Omega1 - Qf1)/N) + N * (sin(beta) * p - Omega1));
f_roll_values = f_roll(N, R);
[max_value, max_index] = max(f_roll_values);
[min_value, min_index] = min(f_roll_values);

index = 0 ;
if abs(max_value) < abs(min_value)
    index = min_index ;
else
    index = max_index ;
end
max_v = f_roll_values(index);
voltage_circ = max_v - N*(sin(beta)*p(index) - Omega1(index));
% fprintf(">>voltage_circ in roll : %f [kV] \n",abs(voltage_circ)/1000);

max_p = (voltage_circ)^2/R/10^6;


if print_result
    fprintf(">> Maximum voltage in roll : %f [kV] \n",abs(max_v/10^3));
    fprintf(">> Maximum power in roll :   %f [MW] \n",max_p);
end

if do_fig
    figure();
    hold on;
    grid on;
    yyaxis left;
    plot(t, e1, 'LineWidth', 2);
    ylabel('Armature voltage [V]    ', 'Interpreter', 'latex');
    yyaxis right;
    plot(t, i1, 'LineWidth', 2);
    ylabel('Electrical current [A]', 'Interpreter', 'latex');
    xlabel('Time [s]', 'Interpreter', 'latex');
    hold off;
    % Création du dossier si nécessaire
    folder_path = 'figures/step_1_roll';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end
    save_path = fullfile(folder_path, 'current_voltage_plot.pdf');  % Chemin complet du fichier de sauvegarde
        
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Sauvegarder la figure en PDF
    print(gcf, save_path, '-dpdf', '-bestfit');
end

end