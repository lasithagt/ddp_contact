function [x,u]= demo_robot
% A demo of iLQG/DDP with car-parking dynamics
clc;
close all

fprintf(['\nA demonstration of the iLQG algorithm '...
'with car parking dynamics.\n'...
'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

global robot_
addpath('./2-DOF-robot/')

dt = 0.01; % Discretization step 

% initialize the robot
L_1xx=0.1; L_1xy=0; L_1xz=0; L_1yy=0.1; L_1yz=0; L_1zz=0.1; l_1x=0.1; l_1y=0.1; l_1z=0; m_1=0.5; fv_1=0; fc_1=0; 
L_2xx=0.1; L_2xy=0; L_2xz=0; L_2yy=0; L_2yz=0; L_2zz=0.1; l_2x=0.1; l_2y=0.1; l_2z=0; m_2=0.1; fv_2=0; fc_2=0;
parms = [L_1xx, L_1xy, L_1xz, L_1yy, L_1yz, L_1zz, l_1x, l_1y, l_1z, m_1, fv_1, fc_1, L_2xx, L_2xy, L_2xz, L_2yy, L_2yz, L_2zz, l_2x, l_2y, l_2z, m_2, fv_2, fc_2];
robot_ = robot_obj([10,10],parms);

T       = 10;                  % horizon
% Initial and final state vectors:
% q = path_states; 
% xf = repmat([pi/6 pi/3 0 0]', 1, T);

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
full_DDP = false;

% set up the optimization problem
DYNCST   = @(x,u,i) robot_dyn_cst(x,u,full_DDP);
x0       = [pi/6 pi/3 0 0]';           % initial state
u0       = 0.05 + zeros(2,T);    % initial controls
Op.lims  = [-.1 .1;             % wheel angle limits (radians)
             -0.1  0.1];            % acceleration limits (m/s^2)
Op.plot = 1;                    % plot the derivatives as well
Op.lims = [];

% x = zeros(4,200);
% u = ones(2,200);
% [f,c] = car_dyn_cst(x,u,full_DDP);

% robot_dynamics(x,u)

% === run the optimization!
[x,u]= iLQG(DYNCST, x0, u0, Op);



function y = robot_dynamics(x,u)

    global robot_
    % === states and controls:
    % x = [x y t v]' = [x; y; car_angle; front_wheel_velocity]
    % u = [w a]'     = [front_wheel_angle; acceleration]
    final = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;     % h = timestep (seconds)

    % states - velocity
    qd1  = x(3,:,:); 
    qd2  = x(4,:,:); 

    dy = [qd1; qd2; robot_.fkdyn(x,u)]   % change in state
    y  = x + dy * dt;                     % new state

   

function c = robot_cost(x, u)
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns

    final = isnan(u(1,:));
    u(:,final)  = 0;

    cu  = 1e-2*[1 .01];         % control cost coefficients

    cf  = [ .1  .1   1  .3];    % final cost coefficients
    pf  = [.01 .01 .01  1]';    % smoothness scales for final cost

    cx  = 1e-3*[1  1];          % running cost coefficients
    px  = [.1 .1]';             % smoothness scales for running cost

    % control cost
    lu    = cu*u.^2;

    % final cost
    if any(final)
       llf      = cf*sabs(x(:,final),pf);
       lf       = double(final);
       lf(final)= llf;
    else
       lf    = 0;
    end

    % running cost
    lx = cx*sabs(x(1:2,:),px);

    % total cost
    c     = lu + lx + lf;

    function y = sabs(x,p)
    % smooth absolute-value function (a.k.a pseudo-Huber)
    y = pp( sqrt(pp(x.^2,p.^2)), -p);


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = robot_dyn_cst(x,u,full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = robot_dynamics(x,u);
    c = robot_cost(x,u);
else
    % state and control indices
    ix = 1:4;
    iu = 5:6;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) robot_dynamics(xu(ix,:),xu(iu,:));
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    N_J = size(J);
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
    xu_cost = @(xu) robot_cost(xu(ix,:),xu(iu,:));
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