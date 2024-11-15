%% AERO0036-1 Spacecraft control project: Design of a spacecraft attitude control system.

clear variables
close all
clc
set(0,'DefaultTextInterpreter','Latex');


do_fig = 0;       % Set to 1 to plot figures
print_result = 1; % Set to 1 to print results

% Reaction wheel parameters
Omega_max  = 7000;                % Reaction wheels have a maximum speed [RPM]
beta       = deg2rad(63.4);       % Reaction wheels are accommodated with a Beta angle [deg]
R          = 50;                  % Reaction wheel motor has an internal resistance [Ohm]
N          = 1;                   % Reaction wheel motor has torque constant [Nm/A]
c          = 10^(-4);             % Reaction wheel bearings & lubricant have a damping factor [Nm/(rad/s)]
rho_w      = 8000;                % Stainless steel of density [kg/m³]

% % Thruster parameters
% Resize the reaction wheel 
%R = 1;                       
%N = 50;



% this inertia is whitout the reaction wheel
Ixx = 10^6;                         % Spacecraft inertia [kg.m²]
Iyy = 10^6;                         % Spacecraft inertia [kg.m²]
Izz = 2 * 10^6;                     % Spacecraft inertia [kg.m²]
time_step = 10^(-3); % Time step for the simulation [s]


fprintf('#################################################### \n')
fprintf('#                 PRELIMINARY CALCULATIONS:        # \n')
fprintf('#################################################### \n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Roll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_roll      = (0:time_step:5) ;       % Time to be able to change its orientation of 90 degrees in roll (+x) [s]
dx_roll     = 2.5 ;                   % First interval time derivate
degree_roll = 90;            % Rotation angle in roll [deg]


[Iw_roll,max_v_roll,max_p_roll] = Compute_roll(Ixx,beta,R,N,c,time_step,Omega_max,t_roll,dx_roll, degree_roll,do_fig,print_result) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pitch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_p = 2.5;           % Time interval to perform the rotation [s]
dx_pitch = 2.5/2;    % First interval time derivative 
degree_p = 30;       % Rotation to perform [deg]



[Iw_pitch,max_v_pitch,max_power_pitch] = Compute_pitch(Omega_max, beta, R, N, c, Iyy, degree_p, t_p, time_step, dx_pitch, do_fig, print_result) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yaw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_las = 4*10^3;
[theta_yaw, w0, o0] = Compute_yaw(T_las,Izz,do_fig,beta,N,R,print_result);


fprintf('The initial angular velocity is %f rad/s and the initial angular position is %f rad\n', w0, o0)


fprintf('#################################################### \n')
fprintf('#                 LQR Controllers                  # \n')
fprintf('#################################################### \n\n')


R = 2; 
% R = 0.101;                      
N = 169.8485;


w.name = 'wheel';
w.beta = beta;
Iw_roll = Iw_roll * 50;
Iw_pitch = Iw_pitch * 50;
w.Iw_r = Iw_roll; 
w.Iw_p = Iw_pitch;
w.R = R;
w.N = N;
w.RPM_max = 7000;
w.e_max = 10^6 ; % [V]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Roll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Roll LQR controller\n')
fprintf('####################\n')

rot.name = 'Roll';
rot.angle = 90;
rot.overshoot = 0.2; 
rot.accuracy = 0.1;
rot.Tf = 8;
rot.t_goal = 5;
rot.speed = 0;
rot.I = Ixx;
motion = 'roll';
A_roll = [0,1;0,sin(beta)/Ixx*(N^2/R+c)*(-2*sin(beta)-Ixx/(Iw_roll*sin(beta)))];
B_roll = [0;sin(beta)/Ixx*N/R];
C_roll = [1,0];
D_roll = 0;


% initialisation 
n = 2;

% Number of wheels:
nb_wheel = 2;

% Call to the function FindQR(), function that find the optimal value for 
% the Q and R matrix in roll.


theta_roll_final = deg2rad(90);
x_0 = [theta_roll_final,0];

% i_r = FindQR_minPower(A_roll, B_roll, C_roll, D_roll, w, rot, n, nb_wheel, 1, 1);
i_r = FindQR_test(A_roll, B_roll, C_roll, D_roll, w, rot, n, nb_wheel, 1, 1, 1);

Q_roll = [i_r^n,0;0,1];
R_roll = 1/i_r;
K_roll = lqr(A_roll,B_roll,Q_roll,R_roll); 
% Forming the closed-loop system
sys = ss(A_roll-B_roll*K_roll,B_roll,C_roll,D_roll);
% Computing the time evolution
[y,t,x] = initial(sys,x_0,time_step);
% Computing the gain and phase margins
[Gm_roll,Pm_roll,~,~] = margin(sys);
fprintf("\nGain margin: %.2f.\n",Gm_roll);
fprintf("Phase margin : %.2f.\n\n",Pm_roll);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pitch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Pitch LQR controller\n')
fprintf('####################\n')

rot.name      = 'Pitch';
rot.angle     = 30;
rot.speed     = 0;
rot.overshoot = 0.05; 
rot.accuracy  = 2/30;
rot.Tf        = 5;
rot.t_goal    = 2.5;

A_pitch = [0,1;0,sin(beta)/Iyy*(N^2/R+c)*(-2*sin(beta)-Iyy/(Iw_pitch*sin(beta)))];
B_pitch = [0;sin(beta)/Iyy*N/R];
C_pitch = [1,0];

i_pitch = FindQR_test(A_pitch, B_pitch, C_pitch, D_roll, w, rot, n, nb_wheel, 1, 1, 1);
Q_pitch = [i_pitch^n,0;0,1];
R_pitch = 1/i_pitch;
K_pitch = lqr(A_pitch,B_pitch,Q_pitch,R_pitch);
% Forming the closed-loop system
sys     = ss(A_pitch-B_pitch*K_pitch,B_pitch,C_pitch,D_roll);
% Computing the time evolution
[y,t,x] = initial(sys,x_0,time_step);
% Computing the gain and phase margins
[Gm_pitch,Pm_pitch,~,~] = margin(sys);
fprintf("\nGain margin: %.2f.\n",Gm_pitch);
fprintf("Phase margin : %.2f.\n\n",Pm_pitch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yaw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Yaw LQR controller\n')
fprintf('####################\n')

rot.name  = 'Yaw';
rot.angle = o0 * 180/pi;
rot.speed = w0;
rot.overshoot = 0;
rot.accuracy = 0.05;
rot.Tf = 20;
rot.t_goal = 5;

A_yaw = [0,1;0,cos(beta)/Izz*(N^2/R+c)*(-4*cos(beta)-Izz/(Iw_pitch*cos(beta)))];
B_yaw = [0;cos(beta)/Izz*N/R];
C_yaw = [1,0];
nb_wheel = 4;


% q_weight_yaw =0.98*10^9;
% Q_yaw = q_weight_yaw*[1 0 ; 0 0];     

% theta_ref_yaw = 0;
% K_yaw = lqr(A_yaw,B_yaw,Q_yaw,1);

n = 4;
i_yaw = FindQR_test(A_yaw, B_yaw, C_yaw, D_roll, w, rot, n, nb_wheel, 1, 1, 1);
Q_yaw = [i_yaw^n,0;0,1];
R_yaw = 1/i_yaw;
K_yaw = lqr(A_yaw,B_yaw,Q_yaw,R_yaw);
% Forming the closed-loop system
sys = ss(A_yaw-B_yaw*K_yaw,B_yaw,C_yaw,D_roll);

% Computing the time evolution
[y,t,x] = initial(sys,x_0,time_step);
% Computing the gain and phase margins
[Gm_yaw,Pm_yaw,~,~] = margin(sys);
fprintf("\nGain margin: %.2f.\n",Gm_yaw);
fprintf("Phase margin : %.2f.\n\n",Pm_yaw);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PID controllers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_r,denom_r] = ss2tf(A_roll,B_roll,C_roll,D_roll);
H_r = tf(num_r,denom_r);
H_r_scaled = H_r * 1e4;
%sisotool(H_r_scaled)


[num_p,denom_p] = ss2tf(A_pitch,B_pitch,C_pitch,D_roll);
H_p = tf(num_p,denom_p);
H_p_scaled = H_p * 1e4;
%sisotool(H_p_scaled)

[num_y,denom_y] = ss2tf(A_yaw,B_yaw,C_yaw,D_roll);
H_y = tf(num_y,denom_y);
H_y_scaled = H_y * 1e4;
%sisotool(H_y_scaled)

