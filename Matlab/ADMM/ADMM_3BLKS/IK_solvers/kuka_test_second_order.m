function [thetalist, fk_current, dexterity, ratio] = kuka_test_second_order()

    close all;
    
    theta_KUKA = zeros(7,1);
    [Slist_K, g_st] = manipulator_POE; % KUKA
    
    S        = Slist_K;
    pose_K   = FKinSpace(g_st, Slist_K, theta_KUKA');
    eomg     = 0.01;
    ev       = 0.001;
    
    rng default
    thetalist0 = [0,0,0,0,0,0,0]' + 0.1 * rand(7,1);
    thetalistd0 = [0,0,0,0,0,0,0]';
    
    
    % cartesian path
    r = 0.06; t = linspace(0, 2 * pi, 100);
    x = r * sin(t); y = r * cos(t); z = 0.8 * ones(1,numel(t));
    R = eye(3); desired_R_euler = rotm2eul(R)';
    
    % admm terms
    rho = [0 0];
    q_bar  = zeros(7, numel(t));
    qd_bar = zeros(7, numel(t));
    
    fk_desired = [x; y; z; repmat(desired_R_euler, 1, length(t))];
    fk_current = zeros(6, length(t));
    thetalist  = zeros(7, length(t));
    thetalistd = zeros(7, length(t));
    terminate  = @(c, i)c < 0.00001 || i > 10;
    j = 0;
    c = cost_trajectory(fk_desired, fk_current);
    
    while ~terminate(c, j)

        for i = 1:numel(t)
            [thetalist(:,i), thetalistd(:,i), success] = IKinSpace_modified_second_order(S, g_st, RpToTrans(R,[x(i),y(i),z(i)]'), thetalist0, thetalistd0, eomg, ev, rho, q_bar(:,i), qd_bar(:,i));
            
            % success
            % thetalist                 = mod(thetalist, 2 * pi);
            
            thetalist0                = thetalist(:,i);
            fk_                       = FKinSpace(g_st, S, thetalist(:,i));
            fk_current(1:3, i)        = fk_(1:3, 4);
            fk_current(4:6, i)        = rotm2eul(fk_(1:3, 1:3)');
        end
        
        % update terminate conditions
        j = j + 1;
        c = cost_trajectory(fk_desired, fk_current);
        plot3(fk_current(1,:), fk_current(2,:), fk_current(3,:)); 
        grid on
    end
    
    figure(2)
    title('Velocity')
    for i = 1:7
       subplot(4,2,i) 
       plot(thetalistd(i,:))
    end
   
    figure(3)
    title('Joint Position')
    for i = 1:7
       subplot(4,2,i) 
       plot(thetalist(i,:))
    end
    
    function c = cost_trajectory(x_desired, x_trajectory)
        c = sum(sum((x_desired - x_trajectory).^2, 1), 2);
    end

    
end