function [thetalist, thetalistd, fk_current] = kuka_second_order_IK(x_des, q0, qd0, rho, q_bar, qd_bar, is_plot)

    % close all;
    
    [Slist_K, g_st] = manipulator_POE; % KUKA
    
    S        = Slist_K;
    eomg     = 0.001;
    ev       = 0.0001;
    
    thetalist0  = q0;
    thetalistd0 = qd0;
    
    % cartesian path
    %     r = 0.06; t = linspace(0, 2 * pi, 100);
    %     x = r * sin(t); y = r * cos(t); z = 0.8 * ones(1,numel(t));
    %     R = eye(3); desired_R_euler = rotm2eul(R)';
    
    % admm terms
    %     rho = [0 0];
    %     q_bar  = zeros(7, numel(t));
    %     qd_bar = zeros(7, numel(t));
    n = size(x_des, 2);
    
    fk_desired = x_des;
    fk_current = zeros(6, n);
    thetalist  = zeros(7, n);
    thetalistd = zeros(7, n);
    terminate  = @(c, i)c < 0.00001 || i > 1;
    j = 0;
    c = cost_trajectory(fk_desired, fk_current);
    
    thetalist(:,1)  = thetalist0;
    thetalistd(:,1) = thetalistd0;
    fk_                 = FKinSpace(g_st, S, thetalist(:,1));
    fk_current(1:3, 1)  = fk_(1:3, 4);
    fk_current(4:6, 1)  = rotm2eul(fk_(1:3, 1:3)');
    
    while ~terminate(c, j)
        for i = 2:n
            
            R     = eul2rotm([x_des(4,i), x_des(5,i), x_des(6,i)]);
            T_des = RpToTrans(R, x_des(1:3,i));
            
            if (i == 2)
                thetalist0     = thetalist(:,1);
                thetalistd0     = thetalistd(:,1);
            end
            
            [thetalist(:,i), thetalistd(:,i), success] = IKinSpace_modified_second_order(S, g_st, T_des, thetalist0, thetalistd0, eomg, ev, rho, q_bar(:,i), qd_bar(:,i));
            
            % success
            % thetalist                 = mod(thetalist, 2 * pi);

            thetalist0                = thetalist(:,i);
            thetalistd0                = thetalistd(:,i);
            
            fk_                       = FKinSpace(g_st, S, thetalist(:,i));
            fk_current(1:3, i)        = fk_(1:3, 4);
            fk_current(4:6, i)        = rotm2eul(fk_(1:3, 1:3)');
        end
        
        % update terminate conditions
        j = j + 1;
        c = cost_trajectory(fk_desired, fk_current);
        if (is_plot)
            figure(1)
            plot3(fk_current(1,:), fk_current(2,:), fk_current(3,:)); 
            axis([-0.08 0.08 -0.08 0.08 0.7 0.95])
        end
        grid on
        
        
    end
    
%     if (is_plot)
%         figure(2)
%         title('Velocity')
%         for i = 1:7
%            subplot(4,2,i) 
%            plot(thetalistd(i, :))
%         end
% 
%         figure(3)
%         title('Joint Position')
%         for i = 1:7
%            subplot(4,2,i) 
%            plot(thetalist(i,:))
%         end
%     end
    
    function c = cost_trajectory(x_desired, x_trajectory)
        c = sum(sum((x_desired - x_trajectory).^2, 1), 2);
    end

    
end