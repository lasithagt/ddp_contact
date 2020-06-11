function [x,u]= demo_robot_7dof_C_dynamics
% A demo of iLQG/DDP with car-parking dynamics
clc;
close all

fprintf(['\nA demonstration of the iLQG algorithm '...
'with car parking dynamics.\n'...
'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

global robot_ Slist M_K T_el dt
addpath('./7-DOF-robot-C-dynamics/')

dt = 0.01; % Discretization step 
t = linspace(0,2*pi,501);
r = 0.1;
xd_x = r*sin(t);
xd_y = r*cos(t);
xd_z = 0.9*ones(1,numel(t));

global xd q_vec q_0
q_vec = [];

thetalist0   = zeros(7,1) + rand(7,1)*0.1;
q_0          = thetalist0;
[Slist, M_K] = manipulator_POE();
T_el         = @(x) [0 -sin(x(1)) cos(x(1))*sin(x(2));0 cos(x(1)) sin(x(1))*sin(x(2));1 0 cos(x(2))];

% % initialize the robot
L_1xx=0.169; L_1xy=0; L_1xz=0; L_1yy=0.1548396; L_1yz=0; L_1zz=0.024052475; l_1x=0; l_1y=-0.03; l_1z=0.12; m_1=4.50275; fv_1=0; fc_1=0; 
L_2xx=0.07361692375; L_2xy=0; L_2xz=0; L_2yy=0.025943256247; L_2yz=0; L_2zz=0.0596744779975; l_2x=0.0003; l_2y=0.059; l_2z=0.042; m_2=4.50275; fv_2=0; fc_2=0;
L_3xx=0.1334; L_3xy=0; L_3xz=0; L_3yy=0.1257; L_3yz=0; L_3zz=0.0127; l_3x=0.0; l_3y=0.03; l_3z=0.13; m_3=3.0; fv_3=0; fc_3=0;
L_4xx=0.1; L_4xy=0; L_4xz=0; L_4yy=0; L_4yz=0; L_4zz=0.1; l_4x=0.0; l_4y=0.067; l_4z=0.034; m_4=2.61155; fv_4=0; fc_4=0;
L_5xx=0.1; L_5xy=0; L_5xz=0; L_5yy=0; L_5yz=0; L_5zz=0.1; l_5x=0.0001; l_5y=0.021; l_5z=0.076; m_5=1.7; fv_5=0; fc_5=0;
L_6xx=0.1; L_6xy=0; L_6xz=0; L_6yy=0; L_6yz=0; L_6zz=0.1; l_6x=0.0; l_6y=0.0006; l_6z=0.0004; m_6=1.8; fv_6=0; fc_6=0;
L_7xx=0.1; L_7xy=0; L_7xz=0; L_7yy=0; L_7yz=0; L_7zz=0.1; l_7x=0.0; l_7y=0.0; l_7z=0.02; m_7=0.35432; fv_7=0; fc_7=0;

parms = [L_1xx, L_1xy, L_1xz, L_1yy, L_1yz, L_1zz, l_1x, l_1y, l_1z, m_1, fv_1, fc_1,...
         L_2xx, L_2xy, L_2xz, L_2yy, L_2yz, L_2zz, l_2x, l_2y, l_2z, m_2, fv_2, fc_2, ...
         L_3xx, L_3xy, L_3xz, L_3yy, L_3yz, L_3zz, l_3x, l_3y, l_3z, m_3, fv_3, fc_3, ...
         L_4xx, L_4xy, L_4xz, L_4yy, L_4yz, L_4zz, l_4x, l_4y, l_4z, m_4, fv_4, fc_4, ...
         L_5xx, L_5xy, L_5xz, L_5yy, L_5yz, L_5zz, l_5x, l_5y, l_5z, m_5, fv_5, fc_5, ...
         L_6xx, L_6xy, L_6xz, L_6yy, L_6yz, L_6zz, l_6x, l_6y, l_6z, m_6, fv_6, fc_6, ...
         L_7xx, L_7xy, L_7xz, L_7yy, L_7yz, L_7zz, l_7x, l_7y, l_7z, m_7, fv_7, fc_7];

robot_ = robot_obj(parms);

T      = 500;                  
% Initial and final state vectors:
% q = path_states; 
% xf = repmat([pi/6 pi/3 0 0]', 1, T);

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
full_DDP = false;
T0         = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1) 1]'];
[q0,s]     = IKinSpace_modified(Slist, M_K, T0, q_0, 0.001, 0.001);
if any(q0>pi)
    for i=1:7
        if (q0(i)>pi)
            q0(i) = -2*pi + q0(i);
        end
    end
end
s

% set up the optimization problem
DYNCST   = @(x,u,i) robot_dyn_cst(x,u,i,full_DDP);
% x0       = [0.1 0.1 0.6 zeros(1,3) zeros(1,6)]';  % initial state
x0       = [q0' 0 0.1 0.9 0 0 0 zeros(1,7) 0 0 0 0 0 0]';
u0       = 0.0 + zeros(6,T);                      % initial controls
Op.lims  = [-.1 .1;                               % wheel angle limits (radians)
             -0.1  0.1];                          % acceleration limits (m/s^2)
Op.plot  = 1;                                     % plot the derivatives as well
Op.lims  = [];
Op.maxIter = 2;

q_d     = repmat(q0, 1, numel(t));
xd      = [q_d;xd_x;xd_y;xd_z;zeros(9,numel(t))];

% robot_dynamics(x,u)

% === run the optimization!
[x,u]= iLQG(DYNCST, x0, u0, Op);



function y = robot_dynamics(x,u)

    global robot_
    % === states and controls:
    % x = [x y t v]' = [x; y; car_angle; front_wheel_velocity]
    % u = [w a]'     = [front_wheel_angle; acceleration]
    final       = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;                      % h = timestep (seconds)

    % states - velocity
    xd  = x(14:end,:,:); 
    
    dyn = robot_.fkdyn(x,u);
    x(1:7,:) = dyn(1:7,:);
    
    
    dy  = [dyn(8:14,:);xd; dyn(15:end,:)];   % change in state
    
    
    
    y   = x + dy * dt;               % new state

   

function c = robot_cost(x, u, i)
    global xd
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns

    final = isnan(u(1,:));
    u(:,final)  = 0;

    cu  = 1e-3 * [1 1 1 1 1 1];                % control cost coefficients

    %     cf  = 10*[1 1 1 .1 .1 .1 .1 .1 .1 .1 .1 .1];           % final cost coefficients
    %     pf  = 10*[1 1 1  1. 1. 1. .1 .1 .1 .1 .1 .1]';      % smoothness scales for final cost
    % 
    %     cx  = 1*1e-2 * [1 1 1 1 1 1 .1 .1 .1 .1 .1 .1];                % running cost coefficients
    %     px  = 1*[.1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1]';                   % smoothness scales for running cost

    cf  = 5e1*[zeros(1,7) 10 10 0 0 0 0 0 0 0 0 0 0];      % final cost coefficients
    pf  = 4*[zeros(1,7) .1 .1 0 0 0 0 0 0 0 0 0 0]';    % smoothness scales for final cost

    cx  = 5e1*[zeros(1,7) 10 10 0 0 0 0 0 0 0 0 0 0];        % running cost coefficients
    px  = 4*[zeros(1,7) .1 .1 0 0 0 0 0 0 0 0 0 0]';           % smoothness scales for running cost

    % control cost
    lu    = cu * u.^2;
    
    x_d = repmat(xd(:,i), 1,size(x,2)/numel(i));
    
    % final cost
    if any(final)
       llf      = cf * sabs(x(:,final)-x_d(:,final),pf);
       lf       = double(final);
       lf(final)= llf;
    else
       lf    = 0;
    end

    % running cost
    
    %     lx = cx*sabs(x(:,:)-[0.1 0.0 0.5 0 0 0 0 0 0 0 0 0]',px);
    lx = cx * sabs(x(:,:)-x_d, px);
    
    % total cost
    c  = lu + lx + lf;

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
    ix = 1:19;
    iu = 20:25;
    
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
        %JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
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