function [thetalist, success] ...
         = IKinSpace_modified(Slist, M, T, thetalist0, eomg, ev)
% *** CHAPTER 6: INVERSE KINEMATICS ***
% Takes Slist: The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       M: The home configuration of the end-effector,
%       T: The desired end-effector configuration Tsd,
%       thetalist0: An initial guess of joint angles that are close to 
%                   satisfying Tsd,
%       eomg: A small positive tolerance on the end-effector orientation 
%             error. The returned joint angles must give an end-effector 
%             orientation error less than eomg,
%       ev: A small positive tolerance on the end-effector linear position 
%           error. The returned joint angles must give an end-effector 
%           position error less than ev.
% Returns thetalist: Joint angles that achieve T within the specified 
%                    tolerances,
%         success: A logical value where TRUE means that the function found
%                  a solution and FALSE means that it ran through the set 
%                  number of maximum iterations without finding a solution
%                  within the tolerances eomg and ev.
% Uses an iterative Newton-Raphson root-finding method.
% The maximum number of iterations before the algorithm is terminated has 
% been hardcoded in as a variable called maxiterations. It is set to 20 at 
% the start of the function, but can be changed if needed.  
% Example Inputs:
% 
% clear; clc;
% Slist = [[0; 0;  1;  4; 0;    0], ...
%        [0; 0;  0;  0; 1;    0], ...
%        [0; 0; -1; -6; 0; -0.1]];
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
% thetalist0 = [1.5; 2.5; 3];
% eomg = 0.01;
% ev = 0.001;
% [thetalist, success] = IKinSpace(Slist, M, T, thetalist0, eomg, ev)
% 
% Output:
% thetalist =
%    1.5707
%    2.9997
%    3.1415
% success =
%     1

alpha = 4;
thetalist = thetalist0;
i = 0;
maxiterations = 50;
Tsb = FKinSpace(M, Slist, thetalist);
Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T));
err = norm(Vs(1: 3)) > eomg || norm(Vs(4: 6)) > ev;
% err_val_n = norm(Vs(4: 6));
% err_val   = err_val_n;

while err && i < maxiterations
    J = JacobianSpace(Slist, thetalist);
    
    W = eye(7);
%     W = diag([5 1 2 1 2 1 2]);

    if rcond(J*J')<0.01
        JINV = W\J'/(J/W*J'+0.1*eye(6));
    else 
        JINV = W\J'/(J/W*J');
    end
    thetalist_ = thetalist + 1*JINV*Vs; 
    
    
    Tsb = FKinSpace(M, Slist, thetalist_);
    Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T));
    
    % line search
    for j=1:numel(alpha)
        theta_null = alpha(j) * (eye(7) - J(1:3,:)'*JINV(:,1:3)')  * null_space(Slist, thetalist)';
        
        Tsb_ = FKinSpace(M, Slist, thetalist_+theta_null);
        Vs_ = Adjoint(Tsb_) * se3ToVec(MatrixLog6(TransInv(Tsb_) * T));
        
        if (norm(Vs_(4:6))>norm(Vs(4:6)))
            % Vs = Vs;
            thetalist = thetalist_;
%             Vs= Vs_;
%             thetalist = thetalist_+theta_null;
        else
            Vs= Vs_;
            thetalist = thetalist_+theta_null;
            break
        end
        thetalist = mod(thetalist,2*pi);
    
    end
   
%     err_val_n = norm(Vs(4: 6));
%     if (err_val_n > err_val)
%         break
%     end
%     err_val = err_val_n;
    norm(Vs);
    err = norm(Vs(1: 3)) > eomg || norm(Vs(4: 6)) > ev;
    i = i + 1;
end
success = ~ err;
end


function dcdq = null_space(S, q_init)

    f  = @(w)det(J_desired(J(w),w)*0.00 - J(w)*J(w)');
    q  = q_init;
    Jq = J(q);
    d  = det(Jq(1:3,:) * Jq(1:3,:)');
    
    q_min = deg2rad([-169 -119 -169 -119 -169 -119 -179]'); q_max = -q_min;
    q_lim_grad = exp((q - q_min)) + exp((q - q_max));
    %     for i = 1:10
    %         d = d + finite_difference(f,q) * ones(7,1)*0.1;
    %     end
    %     cost = @(q)d - det(J(q_init)*J(q_init)');
    
    dcdq = finite_difference(f, q);% + q_lim_grad';
%     dcdq = q_lim_grad';

    function J_v = J(w)
        J_v = JacobianSpace(S,w);
        J_v = J_v(1:3,:);
    end
    function J_d = J_desired(J,w)
        [U,D,V] = svd(J*J');
        J_d     = U*0.05*eye(3)*V';
    end
    function J = finite_difference(fun, q)
        if nargin < 3
            h = 2^-17;
        end
        n       = numel(q);
        fun_0   = fun(q);
        x_d     = repmat(q,1,n) + 1*h*eye(n);
        J       = zeros(1,n);
        for j = 1:n
            J(j) = (fun(x_d(:,j)) - fun_0) / h;
        end
    end
end

