function ret = robot_model_dyn(x, u, u_c)
    % dynamics for a point mass
    ret_qdd = zeros(6, size(x,2));
    ret     = zeros(6, size(x,2));
    
    for i = 1:size(x,2)
        M    = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        I    = [1 0 0;0 1 0;0 0 1];
             
        xdd  = M \ (u(1:3,i) - u_c(1:3,i));
        wdd  = I \ u(4:6,i);                % I is at the base     
        
        ret_qdd(:,i)  = [xdd;wdd];
        ret(:,i)      = ret_qdd(:,i);
    end

end