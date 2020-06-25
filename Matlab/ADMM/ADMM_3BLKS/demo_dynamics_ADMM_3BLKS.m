function [x,u] = demo_dynamics_ADMM_3BLKS
% A demo for Alternating Direction Method of Multiplier implemented on a
% car-parking problem
clc
close all
addpath('./dynamics/', './dynamics/kuka')

global tool
tool = [0 0 0.136]';

full_DDP = false;
% set up the optimization problem
DYNCST          = @(x,u,rhao,x_bar,c_bar,u_bar,thetalist_bar,thetalistd_bar,i) robot_dyn_cst(x, u, i, rhao, x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar, full_DDP);

T       = 100;              % horizon 
% x0      = [0 0.2 -0.1 0 0 0 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
% u0      = -0.1 + zeros(6,T);     % initial controls

Op.lims  = [0 2*pi;             % wheel angle limits (radians)
             -30  30];            % acceleration limits (m/s^2)
Op.plot  = 1;                    % plot the derivatives as well
Op.maxIter = 5;

global x_des
%% desired path to track
t    = linspace(0, 2*pi, T+1);
r    = 0.06;
xd_x = r * cos(t);
xd_y = r * sin(3*t);
xd_z = (0.8) * ones(1,numel(t));
Rd_r = 0 * ones(1, numel(t));
Rd_p = 0 * ones(1, numel(t));
Rd_y = 0 * ones(1, numel(t));
% [xd_x, xd_y, xd_z] = lissajous_curve(t, 1.1608);

xd_f = -0.0 * sin(1*t) - 5.5;

x_des = [xd_x; xd_y; xd_z; Rd_r; Rd_p; Rd_y];



% calculate the center of curvature.
global RC K
[~,RC,K] = curvature([xd_x', xd_y', xd_z']);
RC(1) = RC(2);
RC(end) = RC(end);
K(1,:) = K(2,:);
K(end,:) = K(end-1,:);
K = K';

RC(isnan(RC)) = 1e10;
K(isnan(K))   = 0;


% desired paths (motin profile an force profile)
global xd Slist M_ inertial_params
g   = [0 0 9.81];
% lbr = importrobot('lbr820.urdf'); % 14 kg payload version of the KUKA LBR iiwa, a 7 degree-of-freedom robot manipulator
% lbr.DataFormat = 'column'; % data format for joint positions, joint velocities and joint torques
% lbr.Gravity = g; % gravity vector has to be flipped if using a flipped urdf as acceleration is defined w.r.t base frame not world frame
% inertial_params = KUKA_Inertial_Params(lbr);

params = load('inertial_params_KUKA.mat');
inertial_params = params.inertial_params;

T0          = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1) 1]'];
theta0      = 0 + ones(7,1)*0.1;
theta0(2)   = 0.2;
theta0(4)   = 0.2;
[Slist, M_] = manipulator_POE();
[q0,s]      = IKinSpace_modified_initial(Slist, M_, T0, theta0, 0.00001, 0.00001);

% if any(q0 > pi)
%     for i = 1:7
%         if (q0(i) > pi)
%             q0(i) = -2*pi + q0(i);
%         end
%     end
% end

if (s == 1)
    fprintf("Inverse Solution Found...")
else
    error("Inverse Solution Not Found...")
end

x0       = [q0' zeros(1,7) 0 0 0]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
u0       = -0. + zeros(7, T);                                                     % initial controls

u0 (3,1) = 0;
q_des    = repmat(q0, 1, numel(t));
xd       = [q_des; zeros(7,numel(t)); zeros(2,numel(t)) ;xd_f];


%% For testing IK
% q_bar  = zeros(7, size(x_des,2));
% qd_bar = zeros(7, size(x_des,2));
% 
% [theta0, ~, ~]  = kuka_second_order_IK(x_des, x0(1:7), x0(8:14), [0;0], q_bar, qd_bar, 0);
% 
% u0       = -G_kuka(inertial_params, theta0)'; 

%% === run the optimization! ===
[x,u]= ADMM_DDP_3BLKS(DYNCST, x0, u0, Op);

% animate the resulting trajectory


function y = robot_dynamics(x, u, i)
    global RC K
    % === states and controls:
    final = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;     % h = timestep (seconds)

    RC_d = repmat(RC(i), 1,size(x,2)/numel(i));
    K_d  = repmat(K(:,i), 1,size(x,2)/numel(i));
    
    % states - velocity
    xd  = x(8:14, :, :);
    xdd = fdyn_dynamics_admm_3blk(x, u, RC_d, K_d);
    
    % modify here
    dy  = [xd; xdd(1:end,:)];     % change in state

    y   = x + dy * dt;   % new state



% function c = robot_cost(x, u, i)
%     global xd
%     % cost function for robot problem
%     % sum of 3 terms:
%     % lu: quadratic cost on controls
%     % lf: final cost on distance from target parking configuration
%     % lx: running cost on distance from origin to encourage tight turns
% 
%     final = isnan(u(1,:));
%     u(:,final)  = 0;
% 
%     cu  = 5e-3 * [ones(1,7)];                        % control cost coefficients
% 
% 
%     cf  = 0*5e-1 * [0.0*ones(1,7) 0.0*ones(1,7) 1 1 10];         % final cost coefficients
%     pf  = 0*4e-1 * [0.0*ones(1,7) 0.0*ones(1,7) .01 .01 .1]';    % smoothness scales for final cost
% 
%     cx  = 5e-1 * [0.0*ones(1,7) 0.1*ones(1,7)  0 0 10];        % running cost coefficients
%     px  = 4e-1 * [0.0*ones(1,7) 0.01.*ones(1,7) .0 .00 .1]';           % smoothness scales for running cost
% 
%     cx_b = 1e2 * [10 10 0];
%     px_b = 1e2 * [10 10 0]';
%     
%     % control cost
%     lu  = cu * u.^2;
% 
%     x_d = repmat(xd(:,i), 1,size(x,2)/numel(i));
%     
%     % final cost
%     if any(final)
%        % fk       = fkine(rb, x(1:2,final));
%        % llf      = cf * sabs(fk(:,end),pf);
%        llf      = cf * sabs(x(:,final)-x_d(:,final),pf);
%        lf       = double(final);
%        lf(final)= llf;
%     else
%        lf    = 0;
%     end
% 
%     lx    = cx * sabs(x(:,:) - x_d, px);
%     
%     % total cost
%     c     = lu + lx + lf;


function [c] = admm_robot_cost(x, u, i, rhao, x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar)
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns
    global xd RC
    m = size(u,1);
    n = size(x,1);

    final = isnan(u(1,:));
    u(:,final)  = 0;
    u_bar(:,final) = 0;
    
    RC_expand = repmat(RC(i)', 1, size(x,2)/numel(i));
    cen   = 0.3 * sum(x(15:17,:).^2,1) ./ RC_expand;
    
    cu  = 5e-5*ones(1,7);         % control cost coefficients

    cf  = 5e-1 * [0.0*ones(1,7) 0.00000*ones(1,7) 0.0 0.0 0.05];        % final cost coefficients
    pf  = 4e-1 * [0.0*ones(1,7) 0.01*ones(1,7) .0 .00 .05 ]';    % smoothness scales for final cost

    cx  = 5e-1 * [0.0*ones(1,7) 0.00000*ones(1,7) 0.0 0.0 0.05];           % running cost coefficients
    px  = 4e-1 * [0.0*ones(1,7) 0.01*ones(1,7) .0 .00 .05]';     % smoothness scales for running cost
    cv  = 5e-5*ones(1,7);
    % control cost

    lu    = cu * u.^2 + (rhao(2)/2) * ones(1,m) * (u-u_bar).^2;

    x_d = repmat(xd(:,i), 1, size(x,2)/numel(i));
    
    % final cost
    if any(final)
       llf      = cf(15:17) * (x(15:17,final)-x_d(15:17,final)).^2;
       lf       = final;
       lf(final)= llf;
    else
       lf    = 0;
    end

%     cx(15:17) * sum((x(15:17,:)-x_d(15:17,:)).^2,2);
    % running cost
    lx     = cx(15:17) * (x(15:17,:)-x_d(15:17,:)).^2 + (rhao(1)/2) * ones(1,n)*(x-x_bar).^2 + ...
        (rhao(3)/2)*(cen-c_bar).^2 + (rhao(5)/2) * ones(1,7)*(x(1:7,:)-thetalist_bar).^2 + ...
        (rhao(4)/2) * ones(1,7)*(x(8:14,:)-thetalistd_bar).^2; 
    
    
    % total cost
    c     = lu + lx + lf; 



% utility functions
function y = sabs(x,p)
% smooth absolute-value function (a.k.a pseudo-Huber)
y = pp( sqrt(pp(x.^2,p.^2)), -p);



function [f, c, fx, fu, fxx, fxu, fuu, cx, cu, cxx, cxu, cuu] = robot_dyn_cst(x, u, i, rhao, x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar, full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives
if nargout == 2
    f = robot_dynamics(x,u,i);
    c = admm_robot_cost(x, u, i, rhao, x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar);
else
    % state and control indices
    ix = 1:17;
    iu = 18:24;
    % dynamics first derivatives
    xu_dyn  = @(xu) robot_dynamics(xu(ix,:),xu(iu,:),i);
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    N_J     = size(J);
    
    % dynamics second derivatives
    if full_DDP
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        if length(N_J) <= 2 
            JJ = reshape(JJ,[4 6 N_J(2)]); 
        else 
            JJ = reshape(JJ, [4 6 N_J(2) N_J(3)]); 
        end
        % JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);    
    else
        [fxx,fxu,fuu] = deal([]);
    end    
    
    % cost first derivatives
    xu_cost = @(xu,x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar) admm_robot_cost(xu(ix,:),xu(iu,:),i,rhao, x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar);
    J       = squeeze(finite_difference2(xu_cost, [x; u],x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    % cost second derivatives
    xu_Jcst = @(xu,x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar) squeeze(finite_difference2(xu_cost, xu,x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar));
    JJ      = finite_difference2(xu_Jcst, [x; u],x_bar, c_bar, u_bar, thetalist_bar, thetalistd_bar);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
%     cx      = diag(rhao(5)/2 * ones(1,7)) * (x(1:7,:)-thetalist_bar);
%     cu      = 
%     cuu     =
%     cxx     =
%     cxu     =
%         
%         lx     = cx * sabs(x(:,:)-x_d, px) + (rhao(1)/2) * ones(1,n)*(x-x_bar).^2 + ...
%         (rhao(3)/2)*(cen-c_bar).^2 + (rhao(5)/2) * ones(1,7)*(x(1:7,:)-thetalist_bar).^2 + ...
%         (rhao(4)/2) * ones(1,7)*(x(8:14,:)-thetalistd_bar).^2; 
%     
    [f,c] = deal([]);
    
end


function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);

function J = finite_difference2(fun, x, bar1,bar2,bar3,bar4,bar5,h)
%% simple finite-difference derivatives
% assumes the function fun() is vectorized
if nargin < 8
    h = 2^-17;
end
[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
B1 = repmat(bar1,1,(n+1));
B2 = repmat(bar2,1,(n+1));
B3 = repmat(bar3,1,(n+1));
B4 = repmat(bar4,1,(n+1));
B5 = repmat(bar5,1,(n+1));
% X1       = reshape(X1, n1, K1*(n1+1));
% X2       = reshape(X2, n1, K1*(n1+1));
Y       = fun(X,B1,B2,B3,B4,B5);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);

function c = tt(a,b)
c = bsxfun(@times,a,b);