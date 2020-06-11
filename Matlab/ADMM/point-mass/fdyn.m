function ret = fdyn(x, u)
    % dynamics for a point mass
    ret_qdd = zeros(3, size(x,2));
    ret_fd  = zeros(1, size(x,2));
    ret = zeros(6, size(x,2));
    
    for i = 1:size(x,2)
        M    = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        xdd  = M \ u(:,i);
        K    = 10; % linear
        ret_fd(:,i)   = K * (x(3,i)-0.1);
        ret_qdd(:,i)  = xdd;
        
        ret(:,i) = [ret_qdd(:,i); 0; 0; ret_fd(:,i)];
    end

end