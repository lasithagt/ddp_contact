function [x,u] = demo_vector_mass_kinematics_ADMM
% A demo for Alternating Direction Method of Multiplier implemented on a
% car-parking problem
clc
close all
addpath('./vector-mass/')

full_DDP = false;

% set up the optimization problem
DYNCST                  = @(x,u,i) robot_dyn_cst(x,u,i,full_DDP);
DYNCST_primal           = @(x,u,rhao,x_bar,u_bar,i) admm_robot_cost(x,u,i,rhao,x_bar,u_bar);

T       = 500;              % horizon 
% x0      = [0 0.2 -0.1 0 0 0 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
% u0      = -0.1 + zeros(6,T);     % initial controls

% u0_bar  = .1*randn(2,T);      % initial controls
Op.lims  = [-pi pi;               % wheel angle limits (radians)
             -8  8];            % acceleration limits (m/s^2)
Op.plot = 1;                    % plot the derivatives as well
Op.maxIter = 10;

% prepare the visualization window and graphics callback
t = linspace(0,2*pi,501);
r = 0.1;
xd_x = r*sin(t);
xd_y = r*cos(t);
xd_z = -0.0*ones(1,numel(t));
xd_f = -1*ones(1,numel(t));

% desired paths (motin profile an force profile)
global xd 

T0          = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1) 1]'];
theta0      = 0 + rand(7,1)*0.0001;
[Slist, M_] = manipulator_POE();
[q0,s]      = IKinSpace_modified(Slist, M_, T0, theta0, 0.001, 0.001);

if any(q0>pi)
    for i=1:7
        if (q0(i)>pi)
            q0(i) = -2*pi + q0(i);
        end
    end
end
s

x0      = [q0' 0 0.08 -0.01 0 0 0 0 0 0 0 0 0 0 0 0.1]';   % states = [position_p, position_w,  velocity_p, velocity_w, force]
u0      = -0.1 + zeros(6,T);                               % initial controls

q_d     = repmat(q0, 1, numel(t));
xd      = [q_d;xd_x;xd_y;xd_z;zeros(11,numel(t));xd_f];

% xd   = [xd_x;xd_y;xd_z;zeros(11,numel(t));xd_f];

% === run the optimization!
[x,u]= ADMM_DDP(DYNCST, DYNCST_primal, x0, u0, Op);

% animate the resulting trajectory


function y = robot_dynamics(x, u)

    % === states and controls:
    % x = [x y t v]' = [x; y; car_angle; front_wheel_velocity]
    % u = [w a]'     = [front_wheel_angle; acceleration]
    final = isnan(u(1,:));
    u(:,final)  = 0;

    % constants
    dt  = 0.01;     % h = timestep (seconds)


    % states - velocity
    qd  = x(14:19,:,:);

    xdd = fdyn_kinematics(x,u);
    dy  = [xdd(1:7,:); qd; xdd(8:end,:)];     % change in state
    y   = x + dy * dt;                        % new state


function c = robot_cost(x, u, i)
    global xd
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns

    final = isnan(u(1,:));
    u(:,final)  = 0;

    cu  = 1e-1*[1 1 1 1 1 1];         % control cost coefficients


%     cf  = 5e-1*[10 10 0 0.1 0.1 0.1 1 1 1 0 0 0 1 1 10];          % final cost coefficients
%     pf  = 4e-1*[.1 .1 .0 0.1 0.1 0.1 1 1 1 0 0 0 .01 .01 .1]';    % smoothness scales for final cost
% 
%     cx  = 5e-1*[10 10 0 0.1 0.1 0.1 1 1 1 0.1 0.1 0.1 1 1 10];        % running cost coefficients
%     px  = 4e-1*[.1 .1 0. 0.1 0.1 0.1 1 1 1 0.1 0.1 0.1 .01 .01 .1]';           % smoothness scales for running cost

    cf  = 0*5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0 0 0 0 0 0 1 1 10];      % final cost coefficients
    pf  = 0*4e-1*[0.0*ones(1,7) .1 .1 .0 0.1 0.1 0.1 0 0 0 0 0 0 .01 .01 .1]';    % smoothness scales for final cost

    cx  = 5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0 0 0 0.0 0.0 0.0 1 1 10];        % running cost coefficients
    px  = 4e-1*[0.0*ones(1,7) .1 .1 0. 0.1 0.1 0.1 0 0 0 0.0 0.0 0.0 .01 .01 .1]';           % smoothness scales for running cost

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
    tool = [0 0 0.2]';

   
    R = eul2rotm(x(4:6,:)', 'ZYX');
    v = R .* tool';
    v = v(:,3,:);
    p = x(1:3,:) - reshape(v,3,[],1);
        
 
    % cost for the base movement
    lx_b = cx_b * sabs(p, px_b);
    
    % cost for the ee movement.
    
    lx = cx * sabs(x(:,:)-x_d, px);
    
    % total cost
    c     = lu + lx + lf ;


function c = admm_robot_cost(x, u, i, rhao, x_bar, u_bar)
    % cost function for robot problem
    % sum of 3 terms:
    % lu: quadratic cost on controls
    % lf: final cost on distance from target parking configuration
    % lx: running cost on distance from origin to encourage tight turns
    global xd
    m = size(u,1);
    n = size(x,1);

    final = isnan(u(1,:));
    u(:,final)  = 0;
    u_bar(:,final) = 0;
    
    cu  = 1e-1*[1 1 1 1 1 1];         % control cost coefficients

%     cf  = 0*5e-1*[10 10 0 0.1 0.1 0.1 1 1 1 0 0 0 1 1 10];      % final cost coefficients
%     pf  = 0*4e-1*[.1 .1 .0 0.1 0.1 0.1 1 1 1 0 0 0 .01 .01 .1]';    % smoothness scales for final cost
% 
%     cx  = 5e-1*[10 10 0 0.1 0.1 0.1 1 1 1 0.1 0.1 0.1 1 1 10];        % running cost coefficients
%     px  = 4e-1*[.1 .1 0. 0.1 0.1 0.1 1 1 1 0.1 0.1 0.1 .01 .01 .1]';           % smoothness scales for running cost

    cf  = 0*5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0 0 0 0 0 0 1 1 10];      % final cost coefficients
    pf  = 0*4e-1*[0.0*ones(1,7) .1 .1 .0 0.1 0.1 0.1 0 0 0 0 0 0 .01 .01 .1]';    % smoothness scales for final cost

    cx  = 5e-1*[0.0*ones(1,7) 10 10 0 0.1 0.1 0.1 0 0 0 0.0 0.0 0.0 1 1 10];        % running cost coefficients
    px  = 4e-1*[0.0*ones(1,7) .1 .1 0. 0.1 0.1 0.1 0 0 0 0.0 0.0 0.0 .01 .01 .1]';           % smoothness scales for running cost

    cx_b = 1e1*[10 10 0];
    px_b = 1e1*[10 10 0]';
    
    % control cost
    lu    = cu*u.^2 + (rhao(2)/2)*ones(1,m)*(u-u_bar).^2;
    %sum(u.*((roll/2)*eye(m)*u),1) - roll*sum(u_bar.*u,1) +...
    %(roll/2)*sum(u_bar.*u_bar,1);%(roll/2)*norm(u - u_bar)^2;

    % final cost
    if any(final)
       llf      = cf*sabs(x(:,final),pf);
       lf       = double(final);
       lf(final)= llf;
    else
       lf    = 0;
    end

    % base of the ee
    tool = [0 0 0.2]';

    R = eul2rotm(x(4:6,:)', 'ZYX');
    v = R .* tool';
    v = v(:,3,:);
    p = x(1:3,:) - reshape(v,3,[],1);
        
    x_d = repmat(xd(:,i), 1,size(x,2)/numel(i));
    
    % cost for the base movement
    lx_b = cx_b * sabs(p, px_b);
    
    
    pa    = FK_SE3(Slist, M, x(1:7,:));
    
    % error from forward function
    err_fk = pa(1:3,end) - x_d();
    
    % running cost
    lx = cx*sabs(x(:,:)-x_d,px) + (rhao(1)/2)*ones(1,n)*(x-x_bar).^2; %sum(x.*((roll/2)*eye(n)*x),1) - roll*sum(x_bar.*x,1) + (roll/2)*sum(x_bar.*x_bar,1);%(roll/2)*norm(x - x_bar)^2;
    
    % total cost
    c     = lu + lx + lf ;


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
    ix = 1:22;
    iu = 23:28;
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
function h = car_plot(x,u)

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