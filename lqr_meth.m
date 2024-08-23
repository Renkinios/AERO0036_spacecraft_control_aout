function [i,u] = lqr_meth(time_interval, A, B, C, D, x_0,wheel_inertias,spacecraft_inertias,motion)

% Parameters :
% -----------
% time_interval: array containing the time interval for the simulation
% A: A matrix from the state-space form
% B: B matrix from the state-space form
% C: C matrix from the state-space form
% D: D matrix from the state-space form
% x_0: array containing the initial conditions of the state, i.e.
%[angle,angular velocity] at t=0s
% Wheel_inertias: array containing the inertias of the wheel, i.e. 
%[I_roll, I_pitch]
% Spacecraft_inertias: array containing the inertias of the spacecraft,
% i.e. [I_xx, I_yy, I_zz]
% Motion: String parameter describing the motion considered, i.e.
% "roll", "pitch", or "yaw"
%

    beta = 63.4*pi/180; % [rad]
    time_step = time_interval(2) - time_interval(1); % [s]
    end_time = time_interval(end); % [s]
    max_voltage = 10^5; % [V]
    max_power = 0.75 * 10^6; % [W]
    max_rot_speed_wheel = 7000/60; % [rps]
    i_array = [];
    u_array = [];
    for i=1:300:1e4
        % Forming Q & R matrices
        if strcmp(motion, 'roll')
            m = 6;
        elseif strcmp(motion, 'pitch')
            m = 3;
        elseif strcmp(motion, 'yaw')
            m = 3;
        end
        Q = [i^m 0; 0 1];                                                                                     
        R = 1/i;
        % Computing the optimal gain
        K = lqr(A,B,Q,R);  
        % Forming the closed-loop system
        sys = ss(A-B*K,B,C,D);
        % Computing the time evolution
        [y,t,x] = initial(sys,x_0,time_interval);
        if strcmp(motion, 'roll')
            % Checking for accuracy
            if abs(y(10/time_step + 1)) <= 0.1 * x_0(1)
                % Checking for overshoot
                if or(min(y) < 0 && abs(min(y)) <= 0.2 * x_0(1), min(y)>0)
                    u = -K*x';
                    u = max(abs(u));
                    % Checking for max. voltage
                    
                    if u <= 2*max_voltage
                        max_rot_speed_spacecraft = max(abs(x(:,2)))/(2*pi); 
                        rot_speed_wheel = spacecraft_inertias(1)*max_rot_speed_spacecraft/(2*wheel_inertias(1)*sin(beta));
                        % Checking for max. wheel rotation speed
                        if rot_speed_wheel <= max_rot_speed_wheel
                            % Storing iterating variables
                            i_array = [i_array,i];
                            u_array = [u_array,u];
                        end
                    end 
                end
                continue
            end
        elseif strcmp(motion, 'pitch')
            % Checking for accuracy
            if abs(y(5/time_step + 1)) <= 0.01 * x_0(1)
                % Checking for overshoot
                if or(min(y) < 0 && abs(min(y)) <= 0.05 * x_0(1), min(y)>0)
                    u = -K*x';
                    u = max(abs(u));
                    % Checking for max. voltage
                    if u <= 2*max_voltage 
                        max_rot_speed_spacecraft = max(abs(x(:,2)))/(2*pi); 
                        rot_speed_wheel = spacecraft_inertias(2)*max_rot_speed_spacecraft/(2*wheel_inertias(2)*sin(beta));
                        % Checking for max. wheel rotation speed
                        if rot_speed_wheel <= max_rot_speed_wheel 
                            i_array = [i_array,i];
                            u_array = [u_array,u];
                        end
                    end 
                end
                continue
            end
        elseif strcmp(motion, 'yaw')
           % Checking for accuracy
            if abs(y(5/time_step + 1)) < 0.05 * x_0(1)
                u = -K*x';
                u = max(abs(u));
                % Checking for max. voltage
                if u <= 4*max_voltage
                    max_rot_speed_spacecraft = max(abs(x(:,2)))/(2*pi); 
                    rot_speed_wheel = spacecraft_inertias(3)*max_rot_speed_spacecraft/(2*(wheel_inertias(1)+wheel_inertias(2))*cos(beta));
                    % Checking for max. wheel rotation speed
                    if rot_speed_wheel <= max_rot_speed_wheel
                        % Storing iterating variables
                       i_array = [i_array,i];
                       u_array = [u_array,u];
                    end
                end 
            end
        end
    end
    
    % Checking if at least one configuration has been found
    if isempty(i_array)
        fprintf("No configuration found in %s\n", motion)
        return 
    else
        % Returning the configuration for which u is the lowest, to
        % minimize the required inertias
        [u_min, i_min] = min(u_array);
        i = i_array(i_min);
        u = u_min;
    end
