function [x,u]= demo_robot_Vector_Mass_dynamics
% A demo of iLQG/DDP with vector mass dynamics
clc;
close all
clear all

global DebugLevel;
DebugLevel = 1;

fprintf('Cartesian space Point Mass joint DDP')

addpath('./vector-mass/', './vector-mass/kuka')

% dt = 0.01; % Discretization step 

T      = 500;                  % horizon
% Initial and final state vectors:

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
full_DDP = false;

% set up the optimization problem
DYNCST  = @(x,u,i) robot_dyn_cst(x,u,i,full_DDP);

% u0(3)   = 0;
% Op.lims = [-.1 .1;                  % wheel angle limits (radians)
%            -0.1  0.1];              % acceleration limits (m/s^2)
Op.plot = 1;                        % plot the derivatives as well
Op.lims = [];
Op.maxIter = 50;

t    = linspace(0,2*pi,501);
r    = 0.1;
% q_d = zeros(7,numel(t));
xd_x = r * sin(t);
xd_y = r * cos(t);
xd_z = -0.0 * ones(1,numel(t)) + 0.9;
xd_f = -0.5 * sin(1*t) - 1.5;
[xd_x, xd_y, xd_z] = lissajous_curve(t, 0.9);


global xd q0 Slist M_ lbr inertial_params
g   = [0 0 9.81];
lbr = importrobot('lbr820.urdf'); % 14 kg payload version of the KUKA LBR iiwa, a 7 degree-of-freedom robot manipulator
lbr.DataFormat = 'column'; % data format for joint positions, joint velocities and joint torques
lbr.Gravity = g; % gravity vector has to be flipped if using a flipped urdf as acceleration is defined w.r.t base frame not world frame
params = load('params.mat');
% inertial_params = KUKA_Inertial_Params(lbr);
inertial_params = params.p;

T0          = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1) 1]'];
theta0      = 0 + rand(7,1)*0.1;
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
    fprintf("\nInverse Solution Found...")
else
    warning("\nInverse Solution Not Found...")
end

x0      = [q0' xd_x(1) xd_y(1) xd_z(1) 0 0 0 zeros(1,7) 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
u0      = -0.1 + zeros(13,T);                                               % initial controls

q_d     = repmat(q0, 1, numel(t));
xd      = [q_d;xd_x;xd_y;xd_z;zeros(3,numel(t));zeros(7,numel(t));zeros(6,numel(t));zeros(2,numel(t));xd_f];
% xd   = [xd_x;xd_y;xd_z;zeros(11,numel(t));xd_f];

% robot_dynamics(x,u)

% === run the optimization!
[x,u] = iLQG(DYNCST, x0, u0, Op);

% plot the figure

function y = robot_dynamics(x,u)

    % === states and controls:
    % x = [x y t v]' = [x; y; car_angle; front_wheel_velocity]
    % u = [w a]'     = [front_wheel_angle; acceleration]
    final = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;     % h = timestep (seconds)

    % states - velocity
    xd  = x(14:26,:,:);

    xdd = fdyn_dynamics(x,u);
    
    % modify here
    dy  = [xd; xdd(1:end,:)];     % change in state
    y   = x + dy * dt;   % new state

   

function c = robot_cost(x, u, i)
    global xd 
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns

    final = isnan(u(1,:));
    u(:,final)  = 0;

    cu  = 1e-1*[1 1 1 1 1 1 1*ones(1,7)];         % control cost coefficients

    cf  = 0*5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 1 1 10];      % final cost coefficients
    pf  = 0*4e-1*[0.0*ones(1,7) .1 .1 .0 0.1 0.1 0.1 0.0*ones(1,7) 0 0 0 0 0 0 .01 .01 .1]';    % smoothness scales for final cost

    cx  = 5e-1*[0.0*ones(1,7) 50 50 0 0.1 0.1 0.1 0.0*ones(1,7) 0.0 0.0 0.0 0.0 0.0 0.0 1 1 10];        % running cost coefficients
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

    %     R = eul2rotm(x(4:6,:)', 'ZYX');
    %     v = R .* tool';
    %     v = v(:,3,:);
    %     p = x(1:3,:) - reshape(v,3,[],1);
    
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
    
    err_fk  = pa - x(8:10,:);
    
    err_fk  = sum(err_fk.^2,1);
     
    % cost for the ee movement.
    lx = cx * sabs(x(:,:)-x_d, px);

    % total cost
    c     = lu + lx + 0*lf + 80*err_fk + 20*err_fk_0;

    function y = sabs(x,p)
    % smooth absolute-value function (a.k.a pseudo-Huber)
    y = pp( sqrt(pp(x.^2,p.^2)), -p);


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = robot_dyn_cst(x,u,i,full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = robot_dynamics(x,u);
    c = robot_cost(x,u,i);
else
    % state and control indices
    ix = 1:29;
    iu = 30:42;
    % dynamics first derivatives
    xu_dyn  = @(xu) robot_dynamics(xu(ix,:),xu(iu,:));
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


% ======== graphics functions ========
function h = robot_plot(x,u)

body        = [0.9 2.1 0.3];           % body = [width length curvature]
bodycolor   = 0.5*[1 1 1];
headlights  = [0.25 0.1 .1 body(1)/2]; % headlights [width length curvature x]
lightcolor  = [1 1 0];
wheel       = [0.15 0.4 .06 1.1*body(1) -1.1 .9];  % wheels = [width length curvature x yb yf]
wheelcolor  = 'k';

h = [];

% make wheels
for front = 1:2
   for right = [-1 1]
      h(end+1) = rrect(wheel,wheelcolor)'; %#ok<AGROW>
      if front == 2
         twist(h(end),0,0,u(1))
      end
      twist(h(end),right*wheel(4),wheel(4+front))
   end
end

% make body
h(end+1) = rrect(body,bodycolor);

% make window (hard coded)
h(end+1) = patch([-.8 .8 .7 -.7],.6+.3*[1 1 -1 -1],'w');

% headlights
h(end+1) = rrect(headlights(1:3),lightcolor);
twist(h(end),headlights(4),body(2)-headlights(2))
h(end+1) = rrect(headlights(1:3),lightcolor);
twist(h(end),-headlights(4),body(2)-headlights(2))

% put rear wheels at (0,0)
twist(h,0,-wheel(5))

% align to x-axis
twist(h,0,0,-pi/2)

% make origin (hard coded)
ol = 0.1;
ow = 0.01;
h(end+1) = patch(ol*[-1 1 1 -1],ow*[1 1 -1 -1],'k');
h(end+1) = patch(ow*[1 1 -1 -1],ol*[-1 1 1 -1],'k');

twist(h,x(1),x(2),x(3))

function twist(obj,x,y,theta)
% a planar twist: rotate object by theta, then translate by (x,y)
i = 1i;
if nargin == 3
   theta = 0;
end
for h = obj
   Z = get(h,'xdata') + i*get(h,'ydata');
   Z = Z * exp(i*theta);
   Z = Z + (x + i*y);
   set(h,'xdata',real(Z),'ydata',imag(Z));
end

function h = rrect(wlc, color)
% draw a rounded rectangle (using complex numbers and a kronecker sum :-)

N        = 25; % number of points per corner

width    = wlc(1);
length   = wlc(2);
curve    = wlc(3);

a        = linspace(0,2*pi,4*N);
circle   = curve*exp(1i*a);
width    = width-curve;
length   = length-curve;
rect1    = diag(width*[1 -1 -1 1] + 1i*length *[1 1 -1 -1]);
rectN    = sum(kron(rect1, ones(1,N)), 1) ;
rr       = circle + rectN;
rr       = [rr rr(1)]; % close the curve

h        = patch(real(rr),imag(rr),color);

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);

function c = tt(a,b)
c = bsxfun(@times,a,b);