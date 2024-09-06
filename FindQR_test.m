%% FindQR Function: 

function [i] = FindQR_test(A, B, C, D, w, rot, n, nb_wheel, print, save, plot_figure)
    time_step = 0.01;
    t = 0:time_step:rot.Tf;    % [s], time vector is fully stabilised with high pointing accuracy
    x0 = [rot.angle*pi/180;rot.speed]; % ICs, we consider that the target is 0� and that we 
                               % start from rad2deg(rot.angle) with no initial angular velocity

    u_store   = [];  % [V], total input voltage (for nb_wheel )
    e_store   = [];  % [V], input voltage (for 1 wheel)
    e_sea_min = [];  % [V], input voltage (for 1 wheel)

    pw_store = []; % [RPM], wheel angular velocity
    y_store  = [];  % [rad], angle of the satellite
    t_store  = [];  % [s], time vector 
    i_store  = [];  % [-], iterator vector
    power_store_min = [] ; 
    power_store = [] ;
    voltage_circ_store = [] ;    
    test = [];
    for i=1:50:1e4
        Q       = [i^n 0;                      % The Q matrix is designed to give 
                    0 1];                      % importance to the precision of the movement                                                                    
        R       = 1/i;                         % Relative importance of the input voltage 
        K       = lqr(A,B,Q,R);                % Computation of the K matrix to modify the system 
        sys     = ss(A-B*K,B,C,D);             % New system
        if strcmp(rot.name, 'Yaw')
            K       = lqr(A,B,Q,R);                % Computation of the K matrix to modify the system 
            sys     = ss(A-B*K,B,C,D);             % New system
            x0 = [0; 0] ;
            pertubation = zeros(size(t));      % Entrée de commande par défaut à 0
            pertubation(t <= 0.5) = 4000;   
            % pertubation(t <= 0.5) = rot.angle*pi/180;    

            [y, t, x] = lsim(sys, pertubation, t, x0);
        else
            [y, t, x] = initial(sys, x0, t);
        end
        % [y,t,x] = initial(sys,x0,t);           % Solving the system

        % Computation of the errors
        % y(end)
        % y(rot.t_goal/time_step + 1)
        if strcmp(rot.name, 'Yaw')
            error_acc = (y(rot.t_goal/time_step + 1) <= 0.005 && y(rot.t_goal/time_step + 1) >= - 0.005);
            % fprintf(y(rot.t_goal/time_step + 1), "\n");
            % fprintf("Rotation cible à l'instant %.2f secondes : %.2f\n", rot.t_goal / time_step, y(rot.t_goal/time_step + 1));
            % fprintf("X0 : %.2f\n", x0(1));

            % fprintf("Error acc = %d\n", error_acc);
        else

            if rot.accuracy ~= 0
                error_acc = (y(rot.t_goal/time_step + 1) <= rot.accuracy* x0(1) && y(rot.t_goal/time_step + 1) >= -rot.accuracy* x0(1));
            else 
                error_acc = 1;
            end
        end
            
        if rot.overshoot ~= 0
            % error_ov = (y(rot.t_goal/time_step + 1) <= 0 && y(rot.t_goal/time_step + 1) >= -rot.overshoot* x0(1));
            error_ov = or(min(y) < 0 && abs(min(y)) <= rot.overshoot * x0(1), min(y)>0);
        else 
            error_ov = 1;
        end
        % if strcmp(rot.name, 'Yaw')
        %     fprintf("Error acc = %d\n", error_acc);
        %     fprintf("Error ov = %d\n", error_ov);
        % end
        % If the solution computed is not fairly accurate and do not 
        % overshoot, the configuration is not retained 
        if error_acc && error_ov 
            % if strcmp(rot.name, 'Yaw')
            %     fprintf(">> Configuration found for %s for i = %d.\n",w.name,i);
            % end
            u = -K*x';
            [max_value, max_index] = max(u);
            [min_value, min_index] = min(u);
            index = 0 ;
            if abs(max_value) < abs(min_value)
                e     = min_value/nb_wheel;
                index = min_index;
                % u_store = [u_store, abs(min_value)]; 
            else
                e     = max_value/nb_wheel;
                index = max_index;
                % u_store = [u_store, abs(max_value)]; 
            end
            % If the voltage needed is grater than the maximal one, the
            % configuration is not retained               
            % Computation of the maximal rotation speed            
            % Computation of the maximal rotation speed 
            p = max(abs(x(:,2)))*30/pi; %  --> rad/s to RPM spacecraft speed
            p_n = x(:,2) ; % c'est le vecteur par rapport au temps mais je pense pas que sa soit l'acceleration
            % Depending on the rotation considered, the speed of 
            % rotation of the wheel is calculated
            if strcmp(rot.name, 'Roll') 
                pw    = rot.I*p/(2*w.Iw_r*sin(w.beta));    % of the wheel scalaire acceleration RPM
                omega =  rot.I*p_n/(2*w.Iw_r*sin(w.beta)); % vecteur of the wheel
            elseif strcmp(rot.name, 'Pitch')
                pw    = rot.I*p/(2*w.Iw_p*sin(w.beta));
                omega =  rot.I*p_n/(2*w.Iw_p*sin(w.beta));
            elseif strcmp(rot.name, 'Yaw')
                pw    = rot.I*p/(2*(w.Iw_p + w.Iw_r)*cos(w.beta));
                omega =  rot.I*p_n/(2*(w.Iw_p + w.Iw_r)*cos(w.beta));
            end
            % Computation of the maximal power
            voltage_circ    = e -  w.N * sin(w.beta) * p_n + omega;
            power           = max(voltage_circ)^2/w.R;
            power_store_min = [power_store, power];
            e_sea_min       = [e_sea_min, abs(e)];
            test_power      = e - w.N * sin(w.beta) * p * pi/30 + pw * pi / 30;
            % If the rotation speed of the wheel is grater than the
            % maximal rotation speed, the configuration is not
            % selected.
            %fprintf(">> Input voltage e = %f [V].\n",pw); 
            if abs(e) <= w.e_max
                if pw <= w.RPM_max
                
                    % Store the voltage required for 1 wheel
                    % Store the wheel rotation speed in RPM
                    % Store the solution y for the corresponding voltage
                    % and the corresponding time vector
                    y_store     = [y_store, y]; 
                    t_store     = [t_store, t];
                    % Store the current iteration
                    i_store     = [i_store, i]; 
                    e_store     = [e_store, abs(e)];
                    pw_store    = [pw_store, pw];
                    power_store = [power_store, power];
                    voltage_circ_store = [voltage_circ_store, voltage_circ];
                    test = [test, test_power] ;
                end 
            end 
        end
    end 
    % If e_store is empty that means that no configuration was found. 
    if isempty(e_store) 
        fprintf("No configuration found.\n"); 
        % If u_store is empty that means that it is impossible to achieve 
        % required precision.
        if isempty(e_sea_min)
            fprintf("Impossible to achieve required precision.\n"); 
        % Otherwise the problem comes from the fact that the desired voltage 
        % is too low. We have therefore to increase the desired voltage
        % range.
        else
            [min_volt, Argument] = min(e_sea_min);
            fprintf("Minimal input voltage e = %f [kV].\n",min_volt/10^3);
            fprintf(">> Min power = %f [MW].\n", power_store_min(Argument)/10^6);
        end
    % In this case, one or more configurations are found. The chosen one is 
    % the one that has the maximal input voltage.
    else

        [e, index] = min(e_store);
        % index = (e==max(e));

        y = y_store(:,index)/pi*180;
        y_new = rot.angle*ones(size(y))-y;
        t = t_store(:,index); 
        fprintf(">> Length vecteur y = %d\n",length(i));
        i = i_store(index); 


        pw = pw_store(index);

        power = power_store(index);

        vol_circ = voltage_circ_store(:,index);
        [vol_circ, index_vol_circ] = max(vol_circ);

        test = test(index);

        fprintf(">> Configuration found for %s for i = %d.\n",w.name,i);
        fprintf(">> Input voltage e = %f [kV].\n",e/10^3); 
        % fprintf(">> %s angle after %d [s]: %f [�].\n",rot.name,rot.Tf,rot.angle - y(end)); 
        fprintf(">> Wheel maximum rotation speed: %f [RPM].\n", pw); 
        fprintf(">> Maximum power:    %f [MW].\n", power/10^6);
        fprintf(">> Maximum voltage:  %f [V].\n", vol_circ);
        fprintf(">> Roll rate:        %f [V].\n", test);
        
        % Plot of the corresponding motion
        if plot_figure
            if strcmp(rot.name, 'Yaw')
                figure();
                hold on;
                grid on;
                y_plot = plot(t, y,'b','linewidth',1);
                lgd = y_plot; 

                if rot.accuracy ~= 0
                    ac_h = plot(t, (1 + rot.accuracy) * rad2deg(x0(1))*ones(1, length(t)),'r','linewidth',1);
                    ac_l = plot(t, (1 - rot.accuracy) * rad2deg(x0(1))*ones(1, length(t)),'k','linewidth',1);
                    lgd = [lgd, ac_h, ac_l];
                end

                if rot.overshoot ~= 0
                    ov = plot(t, (1 + rot.overshoot) * rad2deg(x0(1))*ones(1, length(t)),'g','linewidth',1);
                    lgd = [lgd, ov];
                end
            else
                            figure();
                hold on;
                grid on;
                y_plot = plot(t, y_new,'b','linewidth',1);
                lgd = y_plot; 

                if rot.accuracy ~= 0
                    ac_h = plot(t, (1 + rot.accuracy) * rad2deg(x0(1))*ones(1, length(t)),'r','linewidth',1);
                    ac_l = plot(t, (1 - rot.accuracy) * rad2deg(x0(1))*ones(1, length(t)),'k','linewidth',1);
                    lgd = [lgd, ac_h, ac_l];
                end

                if rot.overshoot ~= 0
                    ov = plot(t, (1 + rot.overshoot) * rad2deg(x0(1))*ones(1, length(t)),'g','linewidth',1);
                    lgd = [lgd, ov];
                end
            end

            xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 14);
            ylabel([rot.name, ' angle [$^\circ$]'], 'Interpreter', 'Latex', 'FontSize', 14);

            lgd_labels = {[rot.name, ' angle']};

            if rot.accuracy ~= 0
                lgd_labels = [lgd_labels, 'Accuracy high boundary', 'Accuracy low boundary']; 
            end 
            if rot.overshoot ~= 0
                lgd_labels = [lgd_labels, 'Overshoot boundary'];
            end

            legend(lgd, lgd_labels, 'Location', 'southeast');
            file_name = ['LQR', rot.name, '_angle'];
            if save
                saveas(gcf,file_name,'eps');
            end
        end
    end 
end 