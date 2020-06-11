% This function takes in current q position, new_cartesian point and vel_d
function [q_ret] = IK_SE3(rb, new_cart_pose, guess) 

    poses = new_cart_pose;
    
    p = reshape(poses,4,4,[]);
%     theta_hat = ikine2(rb, p);
%     q_ret     = theta_hat;
%     
    
%     guess = zeros(2,1) + 0.1*rand(2,1);
    
    n = size(p,3);
    for i=1:n
        p_ = twistcoords(logm(p(:,:,i)));
        theta_hat = ikine(rb, p_, guess);
%         theta_hat = ikine2(rb, p_, theta_hat);
        q_ret(:,i) = theta_hat;
        guess = theta_hat;
    end
    
end

%% TODO: optimizes the null space to increase the manipulability and stay
% away from joint limits
function q = redundancy_optimizer(in, q)
% to stay away from joint limits.

%         q_lim_grad = exp((q - input.q_min)) + exp((q - input.q_max))
%         q0 = null(J) * q_lim_grad
end




