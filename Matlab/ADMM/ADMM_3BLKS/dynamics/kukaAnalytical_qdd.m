    %% This function computes the torque output using the 
    function qdd = kukaAnalytical_qdd(params,q_pos,q_vel,tau)
        % States
        qdd = zeros(7, size(q_pos,2));
%         tic
        for i = 1: size(q_pos,2)

            M    = double(M_kuka(params, q_pos'));
            CX2  = double(C_kuka(params, q_pos', q_vel'));
            Grav = -double(G_kuka(params, q_pos'));
            
            % TODO: Add friction model
             qdd(:,i) = (M+0.0*eye(7))\(tau(:, i) - Grav' - CX2'*q_vel); %  INV_kuka(params, q_pos', q_vel', q_accel'); %
        end
%         toc
    end