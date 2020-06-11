function [x,u] = demo_vector_mass_dynamics_ADMM
% A demo for Alternating Direction Method of Multiplier implemented on a
% car-parking problem
clc
close all
addpath('./vector-mass/', './vector-mass/kuka')

full_DDP = false;

global tool
tool = [0 0 0.136]';

% set up the optimization problem
DYNCST                  = @(x,u,i) robot_dyn_cst(x,u,i,full_DDP);
DYNCST_primal           = @(x,u,rhao,x_bar,u_bar,i) admm_robot_cost(x,u,i,rhao,x_bar,u_bar);

T       = 500;              % horizon 
% x0      = [0 0.2 -0.1 0 0 0 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
% u0      = -0.1 + zeros(6,T);     % initial controls

Op.lims  = [-pi pi;               % wheel angle limits (radians)
             -8  8];            % acceleration limits (m/s^2)
Op.plot = 1;                    % plot the derivatives as well
Op.maxIter = 20;

% prepare the visualization window and graphics callback
t = linspace(0,2*pi,501);
r = 0.1;
xd_x = r*sin(t);
xd_y = r*cos(t);
xd_z = 1.1608*ones(1,numel(t));
xd_f = -0.5 * sin(1*t) - 1.5;
[xd_x, xd_y, xd_z] = lissajous_curve(t, 1.1608);

% calculate the center of curvature.
global RC K
[L,RC,K] = curvature([xd_x', xd_y', xd_z']);
RC(1) = RC(2); RC(end) = RC(end);
K(1,:) = K(2,:); K(end,:) = K(end-1,:);
K = K';
% desired paths (motin profile an force profile)
global xd Slist M_ lbr inertial_params
g   = [0 0 9.81];
lbr = importrobot('lbr820.urdf'); % 14 kg payload version of the KUKA LBR iiwa, a 7 degree-of-freedom robot manipulator
lbr.DataFormat = 'column'; % data format for joint positions, joint velocities and joint torques
lbr.Gravity = g; % gravity vector has to be flipped if using a flipped urdf as acceleration is defined w.r.t base frame not world frame
inertial_params = KUKA_Inertial_Params(lbr);

T0          = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1)-tool(3) 1]'];
theta0      = 0 + rand(7,1)*0.1;
theta0(2) = 0.1;
theta0(4) = 0.1;
[Slist, M_] = manipulator_POE();
[q0,s]      = IKinSpace_modified(Slist, M_, T0, theta0, 0.001, 0.001);

if any(q0>pi)
    for i=1:7
        if (q0(i)>pi)
            q0(i) = -2*pi + q0(i);
        end
    end
end

if (s==1)
    fprintf("Inverse Solution Found...")
else
    error("Inverse Solution Not Found...")
end

x0      = [q0' xd_x(1) xd_y(1) xd_z(1) 0 0 0 zeros(1,7) 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
u0      = -0 + zeros(13,T);                                          % initial controls

q_d     = repmat(q0, 1, numel(t));
xd      = [q_d;xd_x;xd_y;xd_z;zeros(3,numel(t));zeros(7,numel(t));zeros(6,numel(t));zeros(2,numel(t));xd_f];
% xd   = [xd_x;xd_y;xd_z;zeros(11,numel(t));xd_f];


% === run the optimization! ===
[x,u]= ADMM_DDP(DYNCST, DYNCST_primal, x0, u0, Op);

% animate the resulting trajectory


function y = robot_dynamics(x, u,i)
    global RC K
    % === states and controls:
    final = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;     % h = timestep (seconds)

    RC_d = repmat(RC(i), 1,size(x,2)/numel(i));
    K_d  = repmat(K(:,i), 1,size(x,2)/numel(i));
    % states - velocity
    xd  = x(14:26,:,:);
    xdd = fdyn_dynamics(x,u,RC_d, K_d);
    
    
    % modify here
    dy  = [xd; xdd(1:end,:)];     % change in state

    y   = x + dy * dt;   % new state


function c = robot_cost(x, u, i)
    global xd tool 
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns

    final = isnan(u(1,:));
    u(:,final)  = 0;

    cu  = 5e-1*[1 1 1 1 1 1 ones(1,7)];         % control cost coefficients


    cf  = 0*5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 1 1 10];      % final cost coefficients
    pf  = 0*4e-1*[0.0*ones(1,7) .1 .1 .0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 .01 .01 .1]';    % smoothness scales for final cost

    cx  = 5e-1*[0.0*ones(1,7) 80 80 0 0.1 0.1 0.1 0.0*ones(1,7) 0.0 0.0 0.0 0.0 0.0 0.0 1 1 10];        % running cost coefficients
    px  = 4e-1*[0.0*ones(1,7) .1 .1 0. 0.1 0.1 0.1 0.*ones(1,7) 0.0 0.0 0.0 0.0 0.0 0.0 .01 .01 .1]';           % smoothness scales for running cost

    cx_b = 1e2*[10 10 0];
    px_b = 1e2*[10 10 0]';
    
    % control cost
    lu  = cu * u.^2;

    x_d = repmat(xd(:,i), 1,size(x,2)/numel(i));
    
    % final cost
    
    if any(final)
       % fk       = fkine(rb, x(1:2,final));
       % llf      = cf * sabs(fk(:,end),pf);
       llf      = cf * sabs(x(:,final)-x_d(:,final),pf);
       lf       = double(final);
       lf(final)= llf;
    else
       lf    = 0;
    end

    % running cost
    %     fk   = fkine(rb, x(1:2,:));
    %     lx   = cx * sabs(fk(1:2,end),px);
    
    % base of the ee

    %     R = eul2rotm(x(4:6,:)', 'ZYX');
    %     v = R .* tool';
    %     v = v(:,3,:);
    %     p = x(1:3,:) - reshape(v,3,[],1);
        
    % cost for the base movement
    %     lx_b = cx_b * sabs(p, px_b);

    % error from forward function
    pa_ = FK_kuka_vec(x(1:7,:));
    
    if (size(pa_,3)>1)
        pa      = permute(pa_(1:3,end,:),[1 3 2]);
    else
        pa      = pa_(1:3,end,:);
    end
    
    rot_eul  = rotm2eul(pa_(1:3,1:3,:),'ZYX');
    rot_eul  = rot_eul';
    temp     = rot_eul.^2;
    err_fk_0 = sum(temp,1); % desired ZYX = [0 0 0];
    
    
    % position at the centroid
    % x_c = x(8:10,:);% - sum(permute(permute(pa_(1:3,1:3,:), [2 1 3]) .* tool,[2 1 3]),2);
    
    % err_fk = sum((pa - x(8:10,:) + sum(permute(permute(pa_(1:3,1:3,:), [2 1 3]) .* tool,[2 1 3]),2)).^2,1);
    % err_fk  = pa - x_c;
    % err_fk  = sum(err_fk.^2,1);
    err_fk  = sum((pa - x(8:10,:) + tool).^2,1);
    
    % cost for the ee movement.
    lx = cx * sabs(x(:,:)-x_d, px);
    
    % total cost
    c     = lu + lx + lf + 1200 * err_fk + 20 * err_fk_0;


function c = admm_robot_cost(x, u, i, rhao, x_bar, u_bar)
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns
    global xd tool
    m = size(u,1);
    n = size(x,1);

    final = isnan(u(1,:));
    u(:,final)  = 0;
    u_bar(:,final) = 0;
    
    cu  = 5e-1*[1 1 1 1 1 1 ones(1,7)];         % control cost coefficients

    cf  = 0*5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 1 1 10];      % final cost coefficients
    pf  = 0*4e-1*[0.0*ones(1,7) .1 .1 .0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 .01 .01 .1]';    % smoothness scales for final cost

    cx  = 5e-1*[0.0*ones(1,7) 80 80 0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0.0 0.0 0.0 1 1 10];        % running cost coefficients
    px  = 4e-1*[0.0*ones(1,7) .1 .1 0. 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0.0 0.0 0.0 .01 .01 .1]';           % smoothness scales for running cost

    % cx_b = 1e1*[10 10 0];
    % px_b = 1e1*[10 10 0]';
    
    % control cost
    lu    = cu*u.^2 + (rhao(2)/2)*ones(1,m)*(u-u_bar).^2;
    % sum(u.*((roll/2)*eye(m)*u),1) - roll*sum(u_bar.*u,1) +...
    % (roll/2)*sum(u_bar.*u_bar,1);%(roll/2)*norm(u - u_bar)^2;

    % final cost
    if any(final)
       llf      = cf*sabs(x(:,final),pf);
       lf       = double(final);
       lf(final)= llf;
    else
       lf    = 0;
    end

    % base of the ee
    %     R = eul2rotm(x(4:6,:)', 'ZYX');
    %     v = R .* tool';
    %     v = v(:,3,:);
    %     p = x(1:3,:) - reshape(v,3,[],1);
        
    x_d = repmat(xd(:,i), 1,size(x,2)/numel(i));
    
    % cost for the base movement
    % lx_b = cx_b * sabs(p, px_b);
    
    
    % error from forward function
    pa_ = FK_kuka_vec(x(1:7,:));
    
    if (size(pa_,3)>1)
        pa      = permute(pa_(1:3,end,:),[1 3 2]);
    else
        pa      = pa_(1:3,end,:);
    end
    
    rot_eul  = rotm2eul(pa_(1:3,1:3,:),'ZYX');
    rot_eul  = rot_eul';
    temp     = rot_eul.^2;
    err_fk_0 = sum(temp,1); % desired ZYX = [0 0 0];
    
    % position at the centroid
    % err_fk = sum((pa -  + sum(permute(permute(pa_(1:3,1:3,:), [2 1 3]) .* tool,[2 1 3]),2)).^2,1);
    
    % err_fk  = pa - x(8:10,:);
    err_fk  = sum((pa - x(8:10,:) + tool).^2,1);
    
    % running cost
    lx = cx*sabs(x(:,:)-x_d,px) + (rhao(1)/2)*ones(1,n)*(x-x_bar).^2; % sum(x.*((roll/2)*eye(n)*x),1) - roll*sum(x_bar.*x,1) + (roll/2)*sum(x_bar.*x_bar,1);%(roll/2)*norm(x - x_bar)^2;
    
    % total cost
    c     = lu + lx + lf + 1200*err_fk + 20*err_fk_0; 


% utility functions
function y = sabs(x,p)
% smooth absolute-value function (a.k.a pseudo-Huber)
y = pp( sqrt(pp(x.^2,p.^2)), -p);



function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = robot_dyn_cst(x,u,i,full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives
if nargout == 2
    f = robot_dynamics(x,u,i);
    c = robot_cost(x,u,i);
else
    % state and control indices
    ix = 1:29;
    iu = 30:42;
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
    xu_cost = @(xu) robot_cost(xu(ix,:),xu(iu,:),i);
    J       = squeeze(finite_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    % cost second derivatives
    xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
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




% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);

function c = tt(a,b)
c = bsxfun(@times,a,b);