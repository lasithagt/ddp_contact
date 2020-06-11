function [thetalist, fk_current, dexterity, ratio] = kuka_test_second_order(theta_KUKA)

    close all;
    
    %% manipulator twists (kuka)
    W1 = [0 0 1]'; W2 = [0 -1 0]';
    w1 = W1; q1 = [0;0;0.2025];
    w2 = W2; q2 = [0;0;0.2025];
    w3 = W1; q3 = [0;0;0.2025];
    w4 = -W2; q4 = [0;0;0.2025+0.42];
    w5 = W1; q5 = [0;0;0.2025+0.42];
    w6 = W2; q6 = [0;0;0.2025+0.42+0.4];
    w7 = W1; q7 = [0;0;0.2025+0.42+0.4+0.126];

    g_st = [1 0 0 0;0 1 0 0;0 0 1 0.2025+0.42+0.4+0.126;0 0 0 1];

    h = 0;
    S1 = ScrewToAxis(q1,w1, h);
    S2 = ScrewToAxis(q2,w2, h);
    S3 = ScrewToAxis(q3,w3, h);
    S4 = ScrewToAxis(q4,w4, h);
    S5 = ScrewToAxis(q5,w5, h);
    S6 = ScrewToAxis(q6,w6, h);
    S7 = ScrewToAxis(q7,w7, h);

    S        = [S1, S2, S3, S4, S5, S6, S7];
    Slist_K  = S;
    pose_K   = FKinSpace(g_st, Slist_K, theta_KUKA);
    eomg     = 0.01;
    ev       = 0.001;
    
    rng default
    thetalist0 = [0,0,0,0,0,0,0]' + 0.1*rand(7,1);
    thetalistd0 = [0,0,0,0,0,0,0]';
    
    % cartesian path
    r = 0.06; t = linspace(0, 2 * pi, 100);
    x = r * sin(t); y = r * cos(t); z = 0.8 * ones(1,numel(t));
    R = eye(3); desired_R_euler = rotm2eul(R)';
    
    fk_desired = [x; y; z; repmat(desired_R_euler, 1, length(t))];
    fk_current = zeros(6, length(t));
    thetalist  = zeros(7, length(t));
    thetalistd = zeros(7, length(t));
    terminate  = @(c, i)c < 0.00001 || i > 10;
    j = 0;
    c = cost_trajectory(fk_desired, fk_current);
    
    while ~terminate(c, j)

        for i = 1:numel(t)
            [thetalist(:,i), thetalistd(:,i), success] = IKinSpace_modified_second_order(S, g_st, RpToTrans(R,[x(i),y(i),z(i)]'), thetalist0, thetalistd0, eomg, ev);
            % success
            thetalist                 = mod(thetalist, 2 * pi);
            J_K_s                     = JacobianSpace(Slist_K, thetalist(:,i)');
            
            % dexterity(i)            = det(J_K_s(1:3,:)*J_K_s(1:3,:)');
            % ratio(i) = max(eig(J_K_s(1:3,:)*J_K_s(1:3,:)')) / min(eig(J_K_s(1:3,:)*J_K_s(1:3,:)'));
            
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
    for i = 1:7
       subplot(4,2,i) 
       plot(thetalistd(i,:))
    end
    
    
    function c = cost_trajectory(x_desired, x_trajectory)
        c = sum(sum((x_desired - x_trajectory).^2, 1), 2);
    end

    
end