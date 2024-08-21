%% AERO0036-1 Spacecraft control project: Design of a spacecraft attitude control system.
do_fig = 0;       % Set to 1 to plot figures
print_result = 1; % Set to 1 to print results

close all
clc
set(0,'DefaultTextInterpreter','Latex');


% Reaction wheel parameters
Omega_max  = 7000;                % Reaction wheels have a maximum speed [RPM]
beta       = deg2rad(63.4);       % Reaction wheels are accommodated with a Beta angle [deg]
R          = 50;                  % Reaction wheel motor has an internal resistance [Ohm]
N          = 1;                   % Reaction wheel motor has torque constant [Nm/A]
c          = 10^(-4);             % Reaction wheel bearings & lubricant have a damping factor [Nm/(rad/s)]
rho_w      = 8000;                % Stainless steel of density [kg/m³]

% % Thruster parameters
% Resize the reaction wheel 
R = 0.5 ;                          % Reaction wheel motor has an internal resistance [Ohm]
N = 145;


% this inertia is whitout the reaction wheel
Ixx = 10^6;                         % Spacecraft inertia [kg.m²]
Iyy = 10^6;                         % Spacecraft inertia [kg.m²]
Izz = 2 * 10^6;                     % Spacecraft inertia [kg.m²]

time_step = 10^(-3); % Time step for the simulation [s]


fprintf('#################################################### \n')
fprintf('#                 PRELIMINARY CALCULATIONS:        # \n')
fprintf('#################################################### \n\n')


[Iw_roll,max_v_roll,max_p_roll] = Compute_roll(Ixx,beta,R,N,c,time_step,Omega_max,do_fig,print_result) ; 


t_pitch = (0:time_step:5); % Time to be able to change its orientation of 30 degrees in pitch (+y) [s]
t_p = 5;           % Time interval to perform the rotation [s]
dx_pitch = 2.5;    % First interval time derivative 
degree_p = 30;     % Rotation to perform [deg]

[Iw_pitch,max_v_pitch,max_power_pitch] = Compute_pitch(Omega_max, beta, R, N, c, Iyy, degree_p, t_p, time_step, dx_pitch, do_fig, print_result) ;