function [x, force, u, cost] = ADMM_CONE(obj,h,DYNCST,x0,u0,COMTask,SwingFootTask,spmode,xf,Op)
%ADMM for three blocks (centroidal_model + whole_body + projection operator)
%the constraints are
%   c        = CoM(q)
%   [m*dc;L] = Ag(q)*dq
%   force    = g_f(q,dq,u)
%   u        = u_bar
%   q        = q_bar
%   force    = force_bar
%---------------------- user-adjustable parameters ------------------------
defaults = {'lims',           [],...            control limits
            'parallel',       true,...          use parallel line-search?
            'Alpha',          10.^linspace(0,-3,11),... backtracking coefficients
            'maxIter',        200,...           maximum iterations        
            'lambda',         1,...             initial value for lambda
            'dlambda',        1,...             initial value for dlambda
            'lambdaFactor',   1.6,...           lambda scaling factor
            'lambdaMax',      1e5,...          lambda maximum value
            'lambdaMin',      1e-5,...          below this value lambda = 0
            'plot',           0,...             0: no;  k>0: every k iters; k<0: every k iters, with derivs window
            'print',          1,...             0: no;  1: final; 2: iter; 3: iter, detailed
            'plotFn',         @(x)0,...         user-defined graphics callback
            'cost',           [],...            initial cost for pre-rolled trajectory            
            };

% --- initial sizes and controls
M   = obj.getMass;
n   = size(x0,1);          % dimension of state vector
m   = size(u0, 1);          % dimension of control vector
N   = size(u0, 2);          % number of state transitions
u   = u0;                   % initial control sequence

% some zero matrices to initialize the trajectory
f_init = zeros(4,N+1);
x_init = zeros(n,N+1);
u_init = zeros(m,N+1);
c_init = zeros(2,N+1);
h_init = zeros(6,N+1);

% --- proccess options
if nargin < 4
    Op = struct();
end
Op  = setOpts(defaults,Op);

verbosity = Op.print;

switch numel(Op.lims)
    case 0
    case 2*m
        Op.lims = sort(Op.lims,2);
    case 2
        Op.lims = ones(m,1)*sort(Op.lims(:))';
    case m
        Op.lims = Op.lims(:)*[-1 1];
    otherwise
        error('limits are of the wrong size')
end

% --- initialize trace data structure
trace = struct('iter',nan,'cost',nan,...
        'dcost',nan);
trace = repmat(trace,[min(Op.maxIter,1e6) 1]);
trace(1).iter = 1;
iter = 1;
% --- initial trajectory
[x,force,un,c0]  = initialtraj(x0,u0,f_init,zeros(1,6),c_init,h_init,f_init,x_init,u_init,f_init,COMTask,SwingFootTask,spmode,xf,DYNCST);
c01 = sum(c0(:));

% [x0b,u0b,~]  = initialtraj(x0,u0_bar,DYNCST);
u = un;

% user plotting
Op.plotFn(x);

% constants, timers, counters
stop        = 0;
dcost       = 0;
print_head  = 6; % print headings every print_head lines
last_head   = print_head;
if verbosity > 0
    fprintf('\n=========== begin ADMM ===========\n');
end

% Initialize dual variebles
roll = [0 1e-2 0 0 0 0]; 
alpha = 1.65;
alphak = 1;
yita = 0.999;
it2 = 10;
% Primal variables
xnew = x;
unew = u;
forcenew = force;
% ufnew = [];

u_bar = zeros(size(u));
x_bar = zeros(size(x));
force_bar = zeros(size(force));
alphak_v = ones(1,it2+1);

u_bar2 = zeros(size(u));
x_bar2 = zeros(size(x)); 
force_bar2 = zeros(size(force));

% Dual variables
c_lambda = zeros(2,size(x,2));
h_lambda = zeros(6,size(x,2));
uf_lambda = zeros(size(force));

x_lambda = zeros(size(x));
u_lambda = zeros(size(u));
force_lambda = zeros(size(force));

x_lambda2 = zeros(size(x));
u_lambda2 = zeros(size(u));
force_lambda2 = zeros(size(force));

% ck = (1/roll)*norm(x_lambda-x_lambda2)^2 + (1/roll)*norm(u_lambda-u_lambda2)^2 + roll*norm(x_bar-x_bar2)^2 + roll*norm(u_bar-u_bar2)^2;
    
res_c = randn(1,it2);
res_h = randn(1,it2);
res_uf = randn(1,it2);
res_clambda = randn(1,it2);
res_mlambda = randn(1,it2);
res_uflambda = randn(1,it2);

res_u = randn(1,it2);
res_x = randn(1,it2);
res_force = randn(1,it2);
res_ulambda = randn(1,it2); 
res_xlambda = randn(1,it2);
res_forcelambda = randn(1,it2);
costcomp = zeros(1, it2);

xcom0 = getCOM(obj,x0(1:6,:));
xc0 = [xcom0(1,1);0;xcom0(2,1);zeros(6,1)];
ufnew3d = load('ForceData/force01.mat');
ufnew3d = ufnew3d.force00;
ufnew3d = [ufnew3d(2,:);zeros(1,N);ufnew3d(1,:)];

% centroidal momentum model (constant case)
[xcnew, ufnew3d] = centroidal_init(obj,h,SwingFootTask,xc0, ufnew3d);

% cnew = xcnew(1:2:3,:);
cnew =load('CoMdata/com01.mat');
cnew = cnew.com;
hnew = [xcnew(7:9,:);M*xcnew(4:6,:)];
ufnew = [ufnew3d([3 1],:);zeros(2,size(ufnew3d,2))];
cnew_old = cnew;
hnew_old = hnew;
ufnew_old = ufnew;

% ufnew3d: [x;y;z]
%% ADMM iteration
for i = 1:it2

    if i<1000
    %% Original simple ADMM 
        %====== centroidal momentum model 
%         [xcnew, ufnew] = centroidal_optimizer(xc0, ufnew, roll4, roll5, roll6);
        [xcnew, ufnew3d] = centroidal_init(obj,h,SwingFootTask,xc0, ufnew3d);
        cnew = xcnew(1:2:3,:);
        hnew = [xcnew(7:9,:);M*xcnew(4:6,:)];
        ufnew = [ufnew3d([3 1],:);zeros(2,size(ufnew3d,2))];
        
        %====== proximal operator to minimize to cost
        [xnew, forcenew, unew, ~] = iLQG_TRK(DYNCST, x0, unew, forcenew, roll, ...
            cnew-c_lambda,hnew-h_lambda,ufnew-uf_lambda, ...
            x_bar-x_lambda,u_bar-u_lambda, force_bar-force_lambda, ...
            COMTask,SwingFootTask,spmode,xf,Op);
        [com, momentum] = wholebodymapping(obj,xnew);
        
        %====== projection operator to satisfy the constraint
        x_bar_old = x_bar;
        u_bar_old = u_bar;
        force_bar_old = force_bar;
        [x_bar, force_bar, u_bar] = proj(xnew+x_lambda, unew+u_lambda, forcenew+force_lambda, 1, Op.lims);

        %====== dual variables update
        c_lambda = c_lambda + com - cnew; %rho1
        h_lambda = h_lambda + momentum - hnew; %rho2
        uf_lambda = uf_lambda + forcenew - ufnew; %rho3
        
        x_lambda = x_lambda + xnew - x_bar; %rho4
        u_lambda = u_lambda + unew - u_bar; %rho5
        force_lambda = force_lambda + forcenew - force_bar; %rho6
        
    %% ADMM with relaxtion 
        %====== proximal operator to minimize to cost
%         [xnew, unew, ~] = iLQG_TRK(DYNCST, DYNCST2, x0, unew, roll, x_bar-x_lambda,u_bar-u_lambda, COMTask,SwingFootTask,spmode,xf,Op);
% 
%         % Relaxtion
%         xnew2 = alpha*xnew + (1-alpha)*x_bar;
%         unew2 = alpha*unew + (1-alpha)*u_bar;
% 
%         %====== project operator to satisfy the constraint
%         x_bar_old = x_bar;
%         u_bar_old = u_bar;
%         [x_bar, u_bar] = proj(xnew2+x_lambda, unew2+u_lambda, Op.lims);
% 
%         %====== dual variables update
%         x_lambda = x_lambda + xnew2 - x_bar;
%         u_lambda = u_lambda + unew2 - u_bar;
        

%%==================================================
    else    
     %% Accelerated ADMM(strongly convex)
%         if i == 100
%             x_bar2 = x_bar;
%             u_bar2 = u_bar;
%             x_lambda2 = x_lambda;
%             u_lambda2 = u_lambda;
%         end
%        
%         %====== proximal operator to minimize to cost
%         [xnew, unew, ~] = iLQG2(DYNCST, DYNCST2, x0, unew, roll, x_bar2-x_lambda2,u_bar2-u_lambda2, Op);
% 
%         %====== project operator to satisfy the constraint
%         x_bar_old = x_bar;
%         u_bar_old = u_bar;
%         [x_bar, u_bar] = proj(xnew+x_lambda2, unew+u_lambda2, Op.lims);
% 
%         %====== dual variables update
%         x_lambda_old = x_lambda;
%         u_lambda_old = u_lambda;
%         x_lambda = x_lambda2 + xnew - x_bar;
%         u_lambda = u_lambda2 + unew - u_bar;
% 
%         %====== Strong convex
%         alphak_v(:,i+1) = (1+sqrt(1+4*alphak_v(:,i)^2))/2;
%         x_bar2 = x_bar + 0.375*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(x_bar-x_bar_old);
%         u_bar2 = u_bar + 0.375*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(u_bar-u_bar_old);
%         x_lambda2 = x_lambda + 0.375*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(x_lambda-x_lambda_old);
%         u_lambda2 = u_lambda + 0.375*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(u_lambda-u_lambda_old);

    %% Accelerated ADMM(weakly convex)
% %         if i == 100
% %             x_bar2 = x_bar;
% %             u_bar2 = u_bar;
% %             x_lambda2 = x_lambda;
% %             u_lambda2 = u_lambda;
% %         end
%         %====== proximal operator to minimize to cost
%         [xnew, unew, ~] = iLQG2(DYNCST, DYNCST2, x0, unew, roll, x_bar2-x_lambda2,u_bar2-u_lambda2, Op);
% 
%     %     % Relaxtion
%     %     xnew2 = alpha*xnew + (1-alpha)*x_bar;
%     %     unew2 = alpha*unew + (1-alpha)*u_bar;
% 
%         %====== project operator to satisfy the constraint
%         x_bar_old = x_bar;
%         u_bar_old = u_bar;
%         [x_bar, u_bar] = proj(xnew+x_lambda2, unew+u_lambda2, Op.lims);
% 
%         %====== dual variables update
%         x_lambda_old = x_lambda;
%         u_lambda_old = u_lambda;
%         x_lambda = x_lambda2 + xnew - x_bar;
%         u_lambda = u_lambda2 + unew - u_bar;
% 
%         %====== Weak convex ()
%         ck_old = ck;
%         ck = (1/roll)*norm(x_lambda-x_lambda2)^2 + (1/roll)*norm(u_lambda-u_lambda2)^2 + roll*norm(x_bar-x_bar2)^2 + roll*norm(u_bar-u_bar2)^2;
% 
%         if ck < yita * ck_old
%             alphak_old = alphak;
%             alphak_v(:,i+1) = (1+sqrt(1+4*alphak_v(:,i)^2))/2;
%             x_bar2 = x_bar +1*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(x_bar-x_bar_old);
%             u_bar2 = u_bar +1*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(u_bar-u_bar_old);
%             x_lambda2 = x_lambda +1*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(x_lambda-x_lambda_old);
%             u_lambda2 = u_lambda +1*((alphak_v(:,i)-1)/alphak_v(:,i+1))*(u_lambda-u_lambda_old);
%         else
%             alphak_v(:,i+1) = 1;
%             x_bar2 = x_bar_old;
%             u_bar2 = u_bar_old;
%             x_lambda2 = x_lambda_old;
%             u_lambda2 = u_lambda_old;
%             ck = (1/yita)*ck_old;
%         end
    end
    %%
    %====== residue
    % primal
    res_c(:,i) = norm(com - cnew);
    res_h(:,i) = norm(momentum - hnew);
    res_uf(:,i) = norm(forcenew - ufnew);
    
    res_u(:,i) = norm(unew - u_bar);
    res_x(:,i) = norm(xnew - x_bar);
    res_force(:,i) = norm(forcenew - force_bar);
    
    %dual
    res_clambda(:,i) = roll(1)*norm(cnew - cnew_old);
    res_mlambda(:,i) = roll(2)*norm(hnew - hnew_old);
    res_uflambda(:,i) = roll(3)*norm(ufnew - ufnew_old);
    
    res_xlambda(:,i) = roll(4)*norm(x_bar - x_bar_old);
    res_ulambda(:,i) = roll(5)*norm(u_bar - u_bar_old);
    res_forcelambda(:,i) = roll(6)*norm(force_bar - force_bar_old);
    
    [~,~,~,cost22]  = initialtraj(x0,unew,f_init,zeros(1,6),c_init,h_init,f_init,x_init,u_init,f_init,...
        COMTask,SwingFootTask,spmode,xf,DYNCST);
    costcomp(:,i) = sum(cost22(:));
    
    %====== varying penalty parameter
    if i > 100
%         if res_c(:,i) > 10*res_clambda(:,i)
%             roll(1) = 2*roll(1);
%             c_lambda = c_lambda/2;
%         elseif res_clambda(:,i) > 10*res_c(:,i)
%             roll(1) = roll(1)/2;
%             c_lambda = c_lambda*2;
%         end
%         
%         if res_h(:,i) > 10*res_mlambda(:,i)
%             roll(2) = 2*roll(2);
%             h_lambda = h_lambda/2;
%         elseif res_mlambda(:,i) > 10*res_h(:,i)
%             roll(2) = roll(2)/2;
%             h_lambda = h_lambda*2;
%         end
%         
%         if res_uf(:,i) > 10*res_uflambda(:,i)
%             roll(3) = 2*roll(3);
%             uf_lambda = uf_lambda/2;
%         elseif res_uflambda(:,i) > 10*res_uf(:,i)
%             roll(3) = roll(3)/2;
%             uf_lambda = uf_lambda*2;
%         end
        
        if res_x(:,i) > 10*res_xlambda(:,i)
            roll(4) = 2*roll(4);
            x_lambda = x_lambda/2;
        elseif res_xlambda(:,i) > 10*res_x(:,i)
            roll(4) = roll(4)/2;
            x_lambda = x_lambda*2;
        end
        
        if res_u(:,i) > 10*res_ulambda(:,i)
            roll(5) = 2*roll(5);
            u_lambda = u_lambda/2;
        elseif res_ulambda(:,i) > 10*res_u(:,i)
            roll(5) = roll(5)/2;
            u_lambda = u_lambda*2;
        end
        
        if res_force(:,i) > 10*res_forcelambda(:,i)
            roll(6) = 2*roll(6);
            force_lambda = force_lambda/2;
        elseif res_forcelambda(:,i) > 10*res_force(:,i)
            roll(6) = roll(6)/2;
            force_lambda = force_lambda*2;
        end
    end
end

figure(15)
ppp = 1:it2+1;
plot(ppp,alphak_v);

%% plot the residue
figure(10)
subplot(1,2,1)
l = 1:it2;
plot(l,res_u,'DisplayName','residue u');
hold on;
plot(l,res_x,'DisplayName','residue x');
% plot(l,res_force,'DisplayName','residue force');
% 
plot(l,res_c,'DisplayName','residue c');
% plot(l,res_h,'DisplayName','residue m');
plot(l,res_uf,'DisplayName','residue uf');

plot(l,res_ulambda,'DisplayName','residue dual u');
plot(l,res_xlambda,'DisplayName','residue dual x');
% plot(l,res_forcelambda,'DisplayName','residue dual force');

title('residue of primal and dual variebles for accelerated ADMM')
xlabel('ADMM iteration')
ylabel('residue')
hold off;

subplot(1,2,2)
jj = 1:it2;
plot(jj,costcomp);
title('cost reduction for accelerated ADMM')
xlabel('ADMM iteration')
ylabel('cost')

figure(11)
ii = 1:N;
subplot(2,2,1)
plot(ii,unew);
hold on 
plot(ii,u_bar);
title('trajectory of control input')
xlabel('time step')
ylabel('torque')
hold off;

subplot(2,2,2)
plot(ii,forcenew(1,:),'DisplayName','new norm');
hold on
plot(ii,forcenew(2,:),'DisplayName','new tang');
plot(ii,ufnew3d(1,:),'DisplayName','standard tang');
plot(ii,ufnew3d(3,:),'DisplayName','standard norm');
title('trajectory of contact force')
xlabel('time step')
ylabel('force')
hold off;

subplot(2,2,3)
aa = 1:N+1;
plot(aa,hnew(2,:));
hold on
plot(aa,momentum(2,:));
title('trajectory of AM')
xlabel('time step')
ylabel('angular momentum')
hold off;

subplot(2,2,4)
aa = 1:N+1;
plot(aa,com);
hold on
plot(aa,cnew);
title('trajectory of CoM')
xlabel('time step')
ylabel('CoM')
hold off;

figure(12)
e1 = (unew(1,:)-u_bar(1,:))./unew(1,:);
e2 = (unew(2,:)-u_bar(2,:))./unew(2,:);
plot(ii,e1);
hold on 
plot(ii,e2);
hold off

[~,~,~,costnew]  = initialtraj(x0,unew,f_init,zeros(1,6),c_init,h_init,f_init,x_init,u_init,f_init,...
    COMTask,SwingFootTask,spmode,xf,DYNCST);

%% ====== STEP 5: accept step (or not)
% print headings
if verbosity >= 1 && last_head == print_head
    last_head = 0;
    fprintf('%-12s','iteration','cost')
    fprintf('\n');
end

% print status
if verbosity >= 1
    fprintf('%-12d%-12.6g%-12.3g\n', ...
        iter, sum(costnew(:)), 0);
    last_head = last_head+1;
end


% accept changes
u              = unew;
x              = xnew;
force          = forcenew;
cost           = costnew;
Op.plotFn(x);
   
    % update trace
%     trace(iter).lambda      = lambda;
%     trace(iter).dlambda     = dlambda;
%     trace(iter).alpha       = alpha;
%     trace(iter).improvement = dcost;
%     trace(iter).cost        = sum(cost(:));
%     trace(iter).reduc_ratio = z;
%     stop = graphics(Op.plot,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace(1:iter),0);

% save lambda/dlambda
if stop
    if verbosity > 0
        fprintf('\nEXIT: Terminated by user\n');
    end
end

if iter == Op.maxIter
    if verbosity > 0
        fprintf('\nEXIT: Maximum iterations reached.\n');
    end
end


function [xnew,forcenew,unew,cnew] = initialtraj(x0,u0,f,roll,c_bar,h_bar,uf_bar, x_bar,u_bar,f_bar,COMTask,SwingFootTask,spmode,xf,DYNCST)
% Generate the initial trajectory 

n        = size(x0,1);
m        = size(u0,1);
N        = size(u0,2);

xnew        = zeros(n,N);
xnew(:,1) = x0(:,1);
unew        = u0;
cnew        = zeros(1,N+1);
forcenew    = zeros(4,N);

for i = 1:N
    [xnew(:,i+1), forcenew(:,i), cnew(:,i)]  = DYNCST(xnew(:,i), unew(:,i),f(:,i),roll,...
        c_bar(:,i),h_bar(:,i),uf_bar(:,i),x_bar(:,i),u_bar(:,i),f_bar(:,i),COMTask(:,i),SwingFootTask(:,i),spmode,xf,i);
end
[~, ~, cnew(:,N+1)] = DYNCST(xnew(:,N+1),nan(m,1),f(:,N+1),roll,...
    c_bar(:,N+1),h_bar(:,N+1),uf_bar(:,N+1),x_bar(:,N+1),u_bar(:,N+1),f_bar(:,N+1),COMTask(:,N+1),SwingFootTask(:,N+1),spmode,xf,i);

function [x2, f2, u2] = proj(xnew, unew, fnew, miu, lims)
% Project operator(control-limit): simply clamp the control output
N = size(unew, 2);
m = size(unew, 1);
n = size(xnew, 1);
u2 = zeros(m,N);
x2 = xnew;
f2 = fnew;

for i =1:N
    
    for j = 4:2:6
        if xnew(j,i+1) > 3.1416
            x2(j,i+1) = 3.1416;
        elseif xnew(j,i+1) < 0
            x2(j,i+1) = 0;
        else
            x2(j,i+1) = xnew(j,i+1);
        end
    end
    
    for k = 1:m
        if unew(k,i) > lims(k,2)
            u2(k,i) = lims(k,2);
        elseif unew(k,i) < lims(k,1)
            u2(k,i) = lims(k,1);
        else
            u2(k,i) = unew(k,i);
        end
    end
    
    for l = 2:2:4
        if fnew(l,i) > miu*abs(fnew(l-1,i))
            f2(l,i) = miu*abs(fnew(l-1,i));
        elseif fnew(l,i) < -miu*abs(fnew(l-1,i))
            f2(l,i) = -miu*abs(fnew(l-1,i));
        else
            f2(l,i) = fnew(l,i);
        end
    end
end


function [com, momentum] = wholebodymapping(obj,xnew)
N = size(xnew,2);
com = zeros(2,N);
momentum = zeros(6,N);
for i = 1:N
    com(:,i) = getCOM(obj,xnew(1:6,i));
    momentum(:,i) = centroidalMomentumMatrix(obj,xnew(1:6,i))*xnew(7:12,i);
end


function [xcnew, ufnew] = centroidal_init(obj, h, SwingFootTask, xc0, uf)
n        = size(xc0,1);
m        = size(uf,1);
N        = size(uf,2);

xcnew        = zeros(n,N);
xcnew(:,1) = xc0;
ufnew        = uf;

for i = 1:N
    xcnew(:,i+1) = centroidal_dynamics(obj,h,SwingFootTask,xcnew(:,i), ufnew(:,i));
end

function y = centroidal_dynamics(obj, h, SwingFootTask, xc, uf)
m    = obj.getMass;
N    = size(xc,2);
c    = xc(1:3,:);
dc   = xc(4:6,:);
L    = xc(7:9,:);
p    = [SwingFootTask(1,1)+0.2;0;0];
ddc  = [uf(1,:)/m; 0; uf(3,:)/m+9.81];
dL   = cross((c-p),uf);

c2   = c + h*dc;
dc2  = dc + h*ddc;
L2   = L + h*dL;
y    = [c2;dc2;L2];

function [x, force, u, L, Vx, Vxx, cost, trace, stop] = iLQG_TRK(DYNCST, ...
    x0, u0, f_new, roll, c_bar, h_bar, uf_bar, x_bar, u_bar, f_bar, ...
    COMTask,SwingFootTask,spmode,xf,Op)
%---------------------- user-adjustable parameters ------------------------
defaults = {'lims',           [],...            control limits
            'parallel',       true,...          use parallel line-search?
            'Alpha',          10.^linspace(0,-3,11),... backtracking coefficients
            'tolFun',         1e-5,...          reduction exit criterion
            'tolGrad',        1e-6,...          gradient exit criterion
            'maxIter',        100,...           maximum iterations            
            'lambda',         1,...             initial value for lambda
            'dlambda',        1,...             initial value for dlambda
            'lambdaFactor',   1.6,...           lambda scaling factor
            'lambdaMax',      1e10,...          lambda maximum value
            'lambdaMin',      1e-6,...          below this value lambda = 0
            'regType',        1,...             regularization type 1: q_uu+lambda*eye(); 2: V_xx+lambda*eye()
            'zMin',           0,...             minimal accepted reduction ratio
            'diffFn',         [],...            user-defined diff for sub-space optimization
            'plot',           0,...             0: no;  k>0: every k iters; k<0: every k iters, with derivs window
            'print',          2,...             0: no;  1: final; 2: iter; 3: iter, detailed
            'plotFn',         @(x)0,...         user-defined graphics callback
            'cost',           [],...            initial cost for pre-rolled trajectory            
            };

% --- initial sizes and controls
n   = size(x0, 1);          % dimension of state vector
nf  = size(f_new, 1);
m   = size(u0, 1);          % dimension of control vector
N   = size(u0, 2);          % number of state transitions
u   = u0;                   % initial control sequence
rho1 = roll(1); %CoM
rho2 = roll(2); %[angular momentum;linear momentum]
rho3 = roll(3); %contact force 
rho4 = roll(4); %joint limit
rho5 = roll(5); %torque limit
rho6 = roll(6); %friction cone 

% --- proccess options
if nargin < 4
    Op = struct();
end
Op  = setOpts(defaults,Op);

Op.print = 2;
verbosity = Op.print;

switch numel(Op.lims)
    case 0
    case 2*m
        Op.lims = sort(Op.lims,2);
    case 2
        Op.lims = ones(m,1)*sort(Op.lims(:))';
    case m
        Op.lims = Op.lims(:)*[-1 1];
    otherwise
        error('limits are of the wrong size')
end

lambda   = Op.lambda;
dlambda  = Op.dlambda;
lims = [];

% --- initialize trace data structure
trace = struct('iter',nan,'lambda',nan,'dlambda',nan,'cost',nan,...
        'alpha',nan,'grad_norm',nan,'improvement',nan,'reduc_ratio',nan,...
        'time_derivs',nan,'time_forward',nan,'time_backward',nan);
trace = repmat(trace,[min(Op.maxIter,1e6) 1]);
trace(1).iter = 1;
trace(1).lambda = lambda;
trace(1).dlambda = dlambda;

% --- initial trajectory
if size(x0,2) == 1
    diverge = true;
    for alpha = Op.Alpha
        [x,force,un,cost]  = forward_pass(x0(:,1),alpha*u,[],[],[],1,...
            DYNCST,lims,[],f_new,roll,c_bar, h_bar, uf_bar, ...
            x_bar, u_bar, f_bar,COMTask,SwingFootTask,spmode,xf);
        % simplistic divergence test
        if all(abs(x(:)) < 1e8)
            u = un;
            diverge = false;
            break
        end
    end
elseif size(x0,2) == N+1 % pre-rolled initial forward pass
    x        = x0;
    diverge  = false;
    if isempty(Op.cost)
        error('pre-rolled initial trajectory requires cost')
    else
        cost     = Op.cost;
    end
else
    error('pre-rolled initial trajectory must be of correct length')
end

trace(1).cost = sum(cost(:));

% user plotting
Op.plotFn(x);

if diverge
    [Vx,Vxx, stop]  = deal(nan);
    L        = zeros(m,n,N);
    cost     = [];
    trace    = trace(1);
    if verbosity > 0
        fprintf('\nEXIT: Initial control sequence caused divergence\n');
    end
    return
end

% constants, timers, counters
flgChange   = 1;
stop        = 0;
dcost       = 0;
z           = 0;
expected    = 0;
print_head  = 6; % print headings every print_head lines
last_head   = print_head;
t_start     = tic;
if verbosity > 0
    fprintf('\n=========== begin iLQG ===========\n');
end
graphics(Op.plot,x,u,cost,zeros(m,n,N),[],[],[],[],[],[],trace,1);
for iter = 1:Op.maxIter
    if stop
        break;
    end
    trace(iter).iter = iter;    
    
    %====== STEP 1: differentiate dynamics and cost along new trajectory
    if flgChange
        t_diff = tic;
        [~,~,~,fx,fu,gx,gu,fxx,fxu,fuu,cx,cu,cf,cff,cxx,cxu,cuu]   = DYNCST(x,[u nan(m,1)],[force nan(nf,1)],...
            roll, c_bar, h_bar, [uf_bar nan(nf,1)], x_bar, [u_bar nan(m,1)], [f_bar nan(nf,1)],...
            COMTask,SwingFootTask,spmode,xf,1:N+1);
%         [cx,cu,cf,cff,cxx,cxu,cuu] = preprocessing(x,u,roll,cx,cu,cxx,cxu,cuu,uf_bar,x_bar,u_bar,f_bar);
        trace(iter).time_derivs = toc(t_diff);
        flgChange   = 0;
    end
    
    %====== STEP 2: backward pass, compute optimal control law and cost-to-go
    backPassDone   = 0;
    while ~backPassDone
        
        t_back   = tic;
        [diverge, Vx, Vxx, l, L, dV] = back_pass(cx,cu,cf,cff,cxx,cxu,cuu,gx,gu,fx,fu,fxx,fxu,fuu,lambda,Op.regType,lims,u);
        trace(iter).time_backward = toc(t_back);
        
        if diverge
            if verbosity > 2
                fprintf('Cholesky failed at timestep %d.\n',diverge);
            end
            dlambda   = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
            lambda    = max(lambda * dlambda, Op.lambdaMin);
            if lambda > Op.lambdaMax
                break;
            end
            continue
        end
        backPassDone      = 1;
    end

    % check for termination due to small gradient
    g_norm         = mean(max(abs(l) ./ (abs(u)+1),[],1));
    trace(iter).grad_norm = g_norm;
    if g_norm < Op.tolGrad && lambda < 1e-5
        dlambda   = min(dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
        lambda    = lambda * dlambda * (lambda > Op.lambdaMin);
        if verbosity > 0
            fprintf('\nSUCCESS: gradient norm < tolGrad\n');
        end
        break;
    end
    
    %====== STEP 3: line-search to find new control sequence, trajectory, cost
    fwdPassDone  = 0;
    if backPassDone
        t_fwd = tic;
        if Op.parallel  % parallel line-search
            [xnew,forcenew,unew,costnew] = forward_pass(x0 ,u, L, x(:,1:N), l,...
                Op.Alpha, DYNCST,lims,Op.diffFn,f_new,roll, ...
                c_bar, h_bar, uf_bar, x_bar, u_bar, f_bar, ...
                COMTask,SwingFootTask,spmode,xf);
            Dcost               = sum(cost(:)) - sum(costnew,2);
            [dcost, w]          = max(Dcost);
            alpha               = Op.Alpha(w);
            expected            = -alpha*(dV(1) + alpha*dV(2));
            if expected > 0
                z = dcost/expected;
            else
                z = sign(dcost);
                warning('non-positive expected reduction: should not occur');
            end
            if (z > Op.zMin)
                fwdPassDone = 1;
                costnew     = costnew(:,:,w);
                xnew        = xnew(:,:,w);
                unew        = unew(:,:,w);
                forcenew    = forcenew(:,:,w);
            end
        else            % serial backtracking line-search
            for alpha = Op.Alpha
                [xnew,forcenew,unew,costnew]   = forward_pass(x0 ,u+l*alpha, L, x(:,1:N),[],1,...
                    DYNCST,lims,Op.diffFn,f_new,roll,...
                    c_bar, h_bar, uf_bar, x_bar, u_bar, f_bar,COMTask,SwingFootTask,spmode,xf);
                dcost    = sum(cost(:)) - sum(costnew(:));
                expected = -alpha*(dV(1) + alpha*dV(2));
                if expected > 0
                    z = dcost/expected;
                else
                    z = sign(dcost);
                    warning('non-positive expected reduction: should not occur');
                end
                if (z > Op.zMin)
                    fwdPassDone = 1;
                    break;
                end
            end
        end
        if ~fwdPassDone
            alpha = nan; % signals failure of forward pass
        end
        trace(iter).time_forward = toc(t_fwd);
    end
    
    %====== STEP 4: accept step (or not), draw graphics, print status
    
    % print headings
    if verbosity > 1 && last_head == print_head
        last_head = 0;
        fprintf('%-12s','iteration','cost','reduction','expected','gradient','log10(lambda)')
        fprintf('\n');
    end
    
    if fwdPassDone
        
        % print status
        if verbosity > 1
            fprintf('%-12d%-12.6g%-12.3g%-12.3g%-12.3g%-12.1f\n', ...
                iter, sum(cost(:)), dcost, expected, g_norm, log10(lambda));
            last_head = last_head+1;
        end
        
        % decrease lambda
        dlambda   = min(dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
        lambda    = lambda * dlambda * (lambda > Op.lambdaMin);
        
        % accept changes
        u              = unew;
        x              = xnew;
        force          = forcenew;
        cost           = costnew;
        flgChange      = 1;
        Op.plotFn(x);
        
        % terminate ?
        if dcost < Op.tolFun
            if verbosity > 0
                fprintf('\nSUCCESS: cost change < tolFun\n');
            end
            break;
        end
        
    else % no cost improvement
        % increase lambda
        dlambda  = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
        lambda   = max(lambda * dlambda, Op.lambdaMin);
        
        % print status
        if verbosity > 1
            fprintf('%-12d%-12s%-12.3g%-12.3g%-12.3g%-12.1f\n', ...
                iter,'NO STEP', dcost, expected, g_norm, log10(lambda));           
            last_head = last_head+1;
        end     
        
        % terminate ?
        if lambda > Op.lambdaMax
            if verbosity > 0
                fprintf('\nEXIT: lambda > lambdaMax\n');
            end
            break;
        end
    end
    % update trace
    trace(iter).lambda      = lambda;
    trace(iter).dlambda     = dlambda;
    trace(iter).alpha       = alpha;
    trace(iter).improvement = dcost;
    trace(iter).cost        = sum(cost(:));
    trace(iter).reduc_ratio = z;
    stop = graphics(Op.plot,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace(1:iter),0);
end

% save lambda/dlambda
trace(iter).lambda      = lambda;
trace(iter).dlambda     = dlambda;

if stop
    if verbosity > 0
        fprintf('\nEXIT: Terminated by user\n');
    end
end

if iter == Op.maxIter
    if verbosity > 0
        fprintf('\nEXIT: Maximum iterations reached.\n');
    end
end


if ~isempty(iter)
    diff_t = [trace(1:iter).time_derivs];
    diff_t = sum(diff_t(~isnan(diff_t)));
    back_t = [trace(1:iter).time_backward];
    back_t = sum(back_t(~isnan(back_t)));
    fwd_t = [trace(1:iter).time_forward];
    fwd_t = sum(fwd_t(~isnan(fwd_t)));
    total_t = toc(t_start);
    if verbosity > 0
        fprintf(['\n'...
            'iterations:   %-3d\n'...
            'final cost:   %-12.7g\n' ...
            'final grad:   %-12.7g\n' ...
            'final lambda: %-12.7e\n' ...
            'time / iter:  %-5.0f ms\n'...
            'total time:   %-5.2f seconds, of which\n'...
            '  derivs:     %-4.1f%%\n'...
            '  back pass:  %-4.1f%%\n'...
            '  fwd pass:   %-4.1f%%\n'...
            '  other:      %-4.1f%% (graphics etc.)\n'...
            '=========== end iLQG ===========\n'],...
            iter,sum(cost(:)),g_norm,lambda,1e3*total_t/iter,total_t,...
            [diff_t, back_t, fwd_t, (total_t-diff_t-back_t-fwd_t)]*100/total_t);
    end
    trace    = trace(~isnan([trace.iter]));
%     timing   = [diff_t back_t fwd_t total_t-diff_t-back_t-fwd_t];
    graphics(Op.plot,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace,2); % draw legend
else
    error('Failure: no iterations completed, something is wrong.')
end


function [cx,cu,cxx,cxu,cuu] = preprocessing(x,u,cx,cu,cxx,cxu,cuu,roll,x_bar,u_bar)
N = size(cx,2);
n = numel(cx)/N;
m = numel(cu)/N;

for i = 1:N-1
    cx(:,i) = cx(:,i) + roll*x(:,i) - roll*x_bar(:,i);
    cu(:,i) = cu(:,i) + roll*u(:,i) - roll*u_bar(:,i);
    cxx(:,:,i) = cxx(:,:,i) + roll*eye(n);
    cuu(:,:,i) = cuu(:,:,i) + roll*eye(m);
    cxu(:,:,i) = cxu(:,:,i);
end
cx(:,N) = cx(:,N) + roll*x(:,N) - roll*x_bar(:,N);
cxx(:,:,N) = cxx(:,:,N) + roll*eye(n);


function [xnew,forcenew,unew,cnew] = forward_pass(x0,u,L,x,du,Alpha,DYNCST,lims,diff,...
    f,roll,c_bar,h_bar,uf_bar,x_bar,u_bar,f_bar,COMTask,SwingFootTask,spmode,xf)
% parallel forward-pass (rollout)
% internally time is on the 3rd dimension, 
% to facillitate vectorized dynamics calls

n        = size(x0,1);
K        = length(Alpha);
K1       = ones(1,K); % useful for expansion
m        = size(u,1);
N        = size(u,2);
nf       = size(f,1);

xnew        = zeros(n,K,N);
xnew(:,:,1) = x0(:,ones(1,K));
unew        = zeros(m,K,N);
cnew        = zeros(1,K,N+1);
forcenew    = zeros(nf,K,N);

for i = 1:N
    unew(:,:,i) = u(:,i*K1);
    
    if ~isempty(du)
        unew(:,:,i) = unew(:,:,i) + du(:,i)*Alpha;
    end    
    
    if ~isempty(L)
        if ~isempty(diff)
            dx = diff(xnew(:,:,i), x(:,i*K1));
        else
            dx          = xnew(:,:,i) - x(:,i*K1);
        end
        unew(:,:,i) = unew(:,:,i) + L(:,:,i)*dx;
    end
    
    if ~isempty(lims)
        unew(:,:,i) = min(lims(:,2*K1), max(lims(:,1*K1), unew(:,:,i)));
    end

    [xnew(:,:,i+1), forcenew(:,:,i), cnew(:,:,i)]  = DYNCST(xnew(:,:,i), unew(:,:,i), f(:,i*K1), roll, ...
        c_bar(:,i*K1),h_bar(:,i*K1),uf_bar(:,i*K1),x_bar(:,i*K1),u_bar(:,i*K1),f_bar(:,i*K1),...
        COMTask(:,i*K1), SwingFootTask(:,i*K1), spmode, xf, i*K1);
%     cnew(:,:,i) = DYNCST2(xnew(:,:,i), unew(:,:,i), roll, x_bar(:,i*K1), u_bar(:,i*K1),COMTask(:,i*K1),SwingFootTask(:,i*K1),spmode,xf,i*K1);
end
[~,~,cnew(:,:,N+1)] = DYNCST(xnew(:,:,N+1), nan(m,K,1), nan(nf,K,1), roll, ...
    c_bar(:,(N+1)*K1),h_bar(:,(N+1)*K1),nan(nf,K,1),x_bar(:,(N+1)*K1),nan(m,K,1),nan(nf,K,1), ...
    COMTask(:,(N+1)*K1),SwingFootTask(:,(N+1)*K1),spmode,xf,i);
% put the time dimension in the columns
xnew = permute(xnew, [1 3 2]);
unew = permute(unew, [1 3 2]);
forcenew = permute(forcenew, [1 3 2]);
cnew = permute(cnew, [1 3 2]);


function [diverge, Vx, Vxx, k, K, dV] = back_pass(cx,cu,cf,cff,cxx,cxu,cuu,gx,gu,fx,fu,fxx,fxu,fuu,lambda,regType,lims,u)
% Perform the Ricatti-Mayne backward pass

% tensor multiplication for DDP terms
vectens = @(a,b) permute(sum(bsxfun(@times,a,b),1), [3 2 1]);

N  = size(cx,2);
n  = numel(cx)/N;
m  = numel(cu)/N;

cx    = reshape(cx,  [n N]);
cu    = reshape(cu,  [m N]);
cxx   = reshape(cxx, [n n N]);
cxu   = reshape(cxu, [n m N]);
cuu   = reshape(cuu, [m m N]);

k     = zeros(m,N-1);
K     = zeros(m,n,N-1);
Vx    = zeros(n,N);
Vxx   = zeros(n,n,N);
dV    = [0 0];

Vx(:,N)     = cx(:,N);
Vxx(:,:,N)  = cxx(:,:,N);

diverge  = 0;
for i = N-1:-1:1
    
    Qu  = cu(:,i)      + gu(:,:,i)'*cf(:,i) + fu(:,:,i)'*Vx(:,i+1);
    Qx  = cx(:,i)      + gx(:,:,i)'*cf(:,i) + fx(:,:,i)'*Vx(:,i+1);
    Qux = cxu(:,:,i)'  + gu(:,:,i)'*cff(:,:,i)*gx(:,:,i) + fu(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
    if ~isempty(fxu)
        fxuVx = vectens(Vx(:,i+1),fxu(:,:,:,i));
        Qux   = Qux + fxuVx;
    end
    
    Quu = cuu(:,:,i)   + gu(:,:,i)'*cff(:,:,i)*gu(:,:,i) + fu(:,:,i)'*Vxx(:,:,i+1)*fu(:,:,i);
    if ~isempty(fuu)
        fuuVx = vectens(Vx(:,i+1),fuu(:,:,:,i));
        Quu   = Quu + fuuVx;
    end
    
    Qxx = cxx(:,:,i)   + gx(:,:,i)'*cff(:,:,i)*gx(:,:,i) + fx(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
    if ~isempty(fxx)
        Qxx = Qxx + vectens(Vx(:,i+1),fxx(:,:,:,i));
    end
    
    Vxx_reg = (Vxx(:,:,i+1) + lambda*eye(n)*(regType == 2));
    
    Qux_reg = cxu(:,:,i)'  + gu(:,:,i)'*cff(:,:,i)*gx(:,:,i) + fu(:,:,i)'*Vxx_reg*fx(:,:,i);
    if ~isempty(fxu)
        Qux_reg = Qux_reg + fxuVx;
    end
    
    QuuF = cuu(:,:,i)  + gu(:,:,i)'*cff(:,:,i)*gu(:,:,i) + fu(:,:,i)'*Vxx_reg*fu(:,:,i) + lambda*eye(m)*(regType == 1);
    
    if ~isempty(fuu)
        QuuF = QuuF + fuuVx;
    end
    
    if nargin < 13 || isempty(lims) || lims(1,1) > lims(1,2)
        % no control limits: Cholesky decomposition, check for non-PD
        [R,d] = chol(QuuF);
        if d ~= 0
            diverge  = i;
            return;
        end
        
        % find control law
        kK = -R\(R'\[Qu Qux_reg]);
        k_i = kK(:,1);
        K_i = kK(:,2:n+1);
        
%     else        % solve Quadratic Program
%         lower = lims(:,1)-u(:,i);
%         upper = lims(:,2)-u(:,i);
%         
%         [k_i,result,R,free] = boxQP(QuuF,Qu,lower,upper,k(:,min(i+1,N-1)));
%         if result < 1
%             diverge  = i;
%             return;
%         end
%         
%         K_i    = zeros(m,n);
%         if any(free)
%             Lfree        = -R\(R'\Qux_reg(free,:));
%             K_i(free,:)   = Lfree;
%         end
        
    end
    
    % update cost-to-go approximation
    dV          = dV + [k_i'*Qu  .5*k_i'*Quu*k_i];
    Vx(:,i)     = Qx  + K_i'*Quu*k_i + K_i'*Qu  + Qux'*k_i;
    Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;
    Vxx(:,:,i)  = .5*(Vxx(:,:,i) + Vxx(:,:,i)');
    
    % save controls/gains
    k(:,i)      = k_i;
    K(:,:,i)    = K_i;
end


function  stop = graphics(figures,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace,init)
stop = 0;

if figures == 0
    return;
end

n  = size(x,1);
N  = size(x,2);
nL = size(L,2);
m  = size(u,1);

cost  = sum(cost,1);
T     = [trace.iter];
T     = T(~isnan(T));
mT    = max(T);

% === first figure
if figures ~= 0  && ( mod(mT,figures) == 0 || init == 2 )
    
    fig1 = findobj(0,'name','iLQG');
    if  isempty(fig1)
        fig1 = figure();
        set(fig1,'NumberTitle','off','Name','iLQG','KeyPressFcn',@Kpress,'user',0,'toolbar','none');
        fprintf('Type ESC in the graphics window to terminate early.\n')
    end
    
    if mT == 1
        set(fig1,'user',0);
    end
    
    set(0,'currentfigure',fig1);
    clf(fig1);
    
    ax1   = subplot(2,2,1);
    set(ax1,'XAxisL','top','YAxisL','right','xlim',[1 N],'xtick',[])
    line(1:N,cost,'linewidth',4,'color',.5*[1 1 1]);
    ax2 = axes('Position',get(ax1,'Position'));
    plot((1:N),x','linewidth',2);
    set(ax2,'xlim',[1 N],'Ygrid','on','YMinorGrid','off','color','none');
    set(ax1,'Position',get(ax2,'Position'));
    double_title(ax1,ax2,'state','running cost')
    
    axL = subplot(2,2,3);
    CO = get(axL,'colororder');
    set(axL,'nextplot','replacechildren','colororder',CO(1:min(n,7),:))
    Lp = reshape(permute(L,[2 1 3]), [nL*m N-1])';
    plot(axL,1:N-1,Lp,'linewidth',1,'color',0.7*[1 1 1]);
    ylim  = get(axL,'Ylim');
    ylim  = [-1 1]*max(abs(ylim));
    set(axL,'XAxisL','top','YAxisL','right','xlim',[1 N],'xtick',[],'ylim',ylim)
    axu = axes('Position',get(axL,'Position'));
    plot(axu,(1:N-1),u(:,1:N-1)','linewidth',2);
    ylim  = get(axu,'Ylim');
    ylim  = [-1 1]*max(abs(ylim));
    set(axu,'xlim',[1 N],'Ygrid','on','YMinorGrid','off','color','none','ylim',ylim);
    set(axL,'Position',get(axu,'Position'));
    double_title(axu,axL,'controls','gains')
    xlabel 'timesteps'
    
    ax1      = subplot(2,2,2);
    set(ax1,'XAxisL','top','YAxisL','right','xlim',[1 mT+eps],'xtick',[])
    hV = line(T,[trace(T).cost],'linewidth',4,'color',.5*[1 1 1]);
    ax2 = axes('Position',get(ax1,'Position'));
    converge = [[trace(T).lambda]' [trace(T).alpha]' [trace(T).grad_norm]' [trace(T).improvement]'];
    hT = semilogy(T,max(0, converge),'.-','linewidth',2,'markersize',10);
    set(ax2,'xlim',[1 mT+eps],'Ygrid','on','YMinorGrid','off','color','none');
    set(ax1,'Position',get(ax2,'Position'));
    double_title(ax1,ax2,'convergence trace','total cost')
    
    subplot(2,2,4);
    plot(T,[trace(T).reduc_ratio]','.-','linewidth',2);
    title 'actual/expected reduction ratio'
    set(gca,'xlim',[0 mT+1],'ylim',[0 2],'Ygrid','on');
    xlabel 'iterations'
    
    set(findobj(fig1,'-property','FontSize'),'FontSize',8)
    stop = get(fig1,'user');
end

if figures < 0  &&  (mod(abs(trace(mT).iter)-1,figures) == 0 || init == 2) && ~isempty(Vx)
    
    fig2 = findobj(0,'name','iLQG - derivatives');
    if  isempty(fig2)
        fig2 = figure();
        set(fig2,'NumberTitle','off','Name','iLQG - derivatives','KeyPressFcn',@Kpress,'user',0);
    end
    
    if length(T) == 1
        set(fig2,'user',0);
    end
    
    set(0,'currentfigure',fig2);
    clf(fig2);
    
    subplot(2,3,1);
    plot(1:N,Vx','linewidth',2);
    set(gca,'xlim',[1 N]);
    title 'V_x'
    grid on;
    
    subplot(2,3,4);
    z = reshape(Vxx,nL^2,N)';
    zd = (1:nL+1:nL^2);
    plot(1:N,z(:,setdiff(1:nL^2,zd)),'color',.5*[1 1 1]);
    hold on;
    plot(1:N,z(:,zd),'linewidth',2);
    hold off
    grid on;
    set(gca,'xlim',[1 N]);
    title 'V_{xx}'
    xlabel 'timesteps'
    
    subplot(2,3,2);
    Nfx     = size(fx,3);
    z = reshape(fx,nL^2,Nfx)';
    zd = (1:n+1:n^2);
    plot(1:Nfx,z(:,setdiff(1:n^2,zd)),'color',.5*[1 1 1]);
    hold on;
    plot(1:Nfx,z,'linewidth',2);
    set(gca,'xlim',[1 Nfx+eps]);
    hold off
    grid on;
    title 'f_{x}'
    
    if numel(fxx) > 0
        fxx = fxx(:,:,:,1:N-1);
        subplot(2,3,5);
        z  = reshape(fxx,[numel(fxx)/(N-1) N-1])';
        plot(1:N-1,z);
        title 'f_{xx}'
        grid on;
        set(gca,'xlim',[1 N-1+eps]);
    end
    
    subplot(2,3,3);
    Nfu     = size(fu,3);
    z = reshape(fu,nL*m,Nfu)';
    plot(1:Nfu,z','linewidth',2);
    set(gca,'xlim',[1 Nfu]);
    title 'f_u'
    grid on;
    
    if numel(fuu) > 0
        subplot(2,3,6);
        fuu = fuu(:,:,:,1:N-1);
        z  = reshape(fuu,[numel(fuu)/(N-1) N-1])';
        plot(1:N-1,z);
        title 'f_{uu}'
        grid on;
        set(gca,'xlim',[1 N-1+eps]);
    end
    
    set(findobj(fig2,'-property','FontSize'),'FontSize',8)
    stop = stop + get(fig2,'user');
end

if init == 1
    figure(fig1);
elseif init == 2
    strings  = {'V','\lambda','\alpha','\partial_uV','\Delta{V}'};
    legend([hV; hT],strings,'Location','Best');
end

drawnow;

function Kpress(src,evnt)
if strcmp(evnt.Key,'escape')
    set(src,'user',1)
end

function double_title(ax1, ax2, title1, title2)

t1 = title(ax1, title1);
set(t1,'units','normalized')
pos1 = get(t1,'position');
t2 = title(ax2, title2);
set(t2,'units','normalized')
pos2 = get(t2,'position');
[pos1(2),pos2(2)] = deal(min(pos1(2),pos2(2)));
pos1(1)  = 0.05;
set(t1,'pos',pos1,'HorizontalAlignment','left')
pos2(1)  = 1-0.05;
set(t2,'pos',pos2,'HorizontalAlignment','right')

% setOpts - a utility function for setting default parameters
% ===============
% defaults  - either a cell array or a structure of field/default-value pairs.
% options   - either a cell array or a structure of values which override the defaults.
% params    - structure containing the union of fields in both inputs.
function params = setOpts(defaults,options)

if nargin==1 || isempty(options)
    user_fields  = [];
else
    if isstruct(options)
        user_fields   = fieldnames(options);
    else
        user_fields = options(1:2:end);
        options     = struct(options{:});
    end
end

if isstruct(defaults)
    params   = defaults;
else
    params   = struct(defaults{:});
end

for k = 1:length(user_fields)
    params.(user_fields{k}) = options.(user_fields{k});
end
