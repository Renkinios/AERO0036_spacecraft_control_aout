function [theta_yaw, w0, o0]= Compute_yaw(T_las,Izz,do_fig,beta,N)

Iw = 500; % à changer quand on sait 


t_las = 0.5; 
r_dot_las = T_las / Izz;
r_las = r_dot_las * t_las; 
theta_yaw = 0.5 * r_dot_las * t_las^2 ; 

H_zz_stop = Izz * r_las; 

t_yaw = 5;
r_dot = (2*(theta_yaw/2))/((t_yaw/2)^2);
r = r_dot*(t_yaw/2); 
H_zz = Izz * r;
T_zz = Izz * r_dot;

w0 = r_las;
o0 = theta_yaw;

fprintf('The initial angular velocity is %f rad/s and the initial angular position is %f rad\n', w0, o0) ;   

f = @(x) [
    0.5*x(3)*x(4)^2 + w0*x(4) + o0 - x(1);
    x(3)*x(4) + w0 - x(2);
    0.5*x(5)*x(6)^2 + x(2)*x(6) + x(1); 
    x(5)*x(6) + x(2);
    x(4) - x(6); % We choose to impose the maneuver time equality
    5 - x(4) - x(6)
];

% Initial guess
x0 = [0, 0, 0, 0, 0, 0];

% Solve of the system of equations 
x = fsolve(f, x0);

% Solution
o1 = x(1);
w1 = x(2);
a1 = x(3); 
t1 = x(4);
a2 = x(5); 
t2 = x(6);

t_yaw_las = 0 : 0.01:t_las;
t_yaw_recover1 = 0.01:0.01:t1;
t_yaw_recover2 = 0.01:0.01:t2;

theta_yaw_las = 0.5*r_dot_las*t_yaw_las.^2; 
w_yaw_las = r_dot_las*t_yaw_las;

theta_yaw_recover1 = 0.5*a1*t_yaw_recover1.^2 + w0*t_yaw_recover1 + w0;
w_yaw_recover1 = a1*t_yaw_recover1 + w0;

theta_yaw_recover2 = 0.5*a2*t_yaw_recover2.^2 + w_yaw_recover1(end)*t_yaw_recover2+theta_yaw_recover1(end);
w_yaw_recover2 = a2*t_yaw_recover2+w_yaw_recover1(end);

t_yaw_plot = 0:0.01:(t_yaw+t_las);
thetha_yaw = [theta_yaw_las theta_yaw_recover1 theta_yaw_recover2];
w_yaw = [w_yaw_las w_yaw_recover1 w_yaw_recover2];
H_yaw = Izz*w_yaw; % [(kg.m^2.rad)/s] Angular momentum
                             
T_yaw = [ones(1,length(t_yaw_las))*r_dot_las*Izz ones(1,length(t_yaw_recover1))*a1*Izz ones(1,length(t_yaw_recover2))*a2*Izz]; % [N.m] Torque

r_dot_wheel1 = a1*Izz/(4*Iw*cosd(beta));    
r_dot_wheel2 = a2*Izz/(4*Iw*cosd(beta));
T_yy_wheel_plot = [ones(1,length(t_yaw_las))*0 ones(1,length(t_yaw_recover1))*r_dot_wheel1*Iw ones(1,length(t_yaw_recover2))*(r_dot_wheel2)*Iw];
i_yaw = T_yy_wheel_plot/N;  % [A] Current in the wheel

%e_w = (i_yaw)^2*R;

if do_fig
    figure
    hold on
    plot(t_yaw_plot, H_yaw, 'linewidth', 2)
    xlabel('Time [s]','interpreter','latex')
    ylabel('Angular momentum [Nm $\cdot$ s]','interpreter','latex')
    xlim([0 5.5])  % Définir la limite de l'axe des x

    grid on
    box on
    
    % Définir la taille de la figure (en pouces)
    width = 6;  % Largeur en pouces
    height = 4; % Hauteur en pouces
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, width, height]);

     % Création du dossier si nécessaire
    folder_path = 'figures/step_1_yaw_';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end

    % Définir le chemin complet pour le fichier PDF
    save_path = fullfile(folder_path, 'Momentum.pdf');  % Chemin complet du fichier de sauvegarde
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [width, height]);
    
    % Sauvegarder la figure en PDF avec un ajustement optimal
    print(gcf, save_path, '-dpdf', '-r0');
end 



% --------------- Torque ------------------------------------------------ %
if do_fig
    figure
    hold on
    stairs(t_yaw_plot, T_yaw, 'linewidth', 2)
    xlabel('Time [s]','interpreter','latex')
    ylabel('Torque [Nm]','interpreter','latex')
    xlim([0 5.5])  % Définir la limite de l'axe des x
    grid on
    box on
    
    % Définir la taille de la figure (en pouces)
    width = 6;  % Largeur en pouces
    height = 4; % Hauteur en pouces
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, width, height]);

     % Création du dossier si nécessaire
    folder_path = 'figures/step_1_yaw_';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end

    % Définir le chemin complet pour le fichier PDF
    save_path = fullfile(folder_path, 'Torque.pdf');  % Chemin complet du fichier de sauvegarde
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [width, height]);
    
    % Sauvegarder la figure en PDF avec un ajustement optimal
    print(gcf, save_path, '-dpdf', '-r0');
    
end
% --------------- Current ----------------------------------------------- %    
if do_fig
    figure
    hold on
    stairs(t_yaw_plot, i_yaw, 'linewidth', 2)
    xlabel('Time [s]','interpreter','latex')
    ylabel('Current [A]','interpreter','latex')
    
    xlim([0 5.5])  % Définir la limite de l'axe des x
    grid on
    box on
    
    % Définir la taille de la figure (en pouces)
    width = 6;  % Largeur en pouces
    height = 4; % Hauteur en pouces
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, width, height]);

     % Création du dossier si nécessaire
    folder_path = 'figures/step_1_yaw_';  % Chemin du dossier
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);  % Créer le dossier s'il n'existe pas
    end

    % Définir le chemin complet pour le fichier PDF
    save_path = fullfile(folder_path, 'Current.pdf');  % Chemin complet du fichier de sauvegarde
    
    % Ajuster les propriétés de la figure et du papier pour une meilleure sauvegarde
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [width, height]);
    
    % Sauvegarder la figure en PDF avec un ajustement optimal
    print(gcf, save_path, '-dpdf', '-r0');
    
end

end
