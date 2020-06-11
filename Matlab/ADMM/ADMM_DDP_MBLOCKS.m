function [x, u, cost] = ADMM_DDP_MBLOCKS(DYNCST, DYNCST_CONTACT, DYNCST_primal,DYNCYT_primal_CONTACT, x0, x0_c, u0, u0_c, Op)
%---------------------- user-adjustable parameters ------------------------

% --- initial sizes and controls
n   = size(x0,1);           % dimension of state vector
n_c = size(x0_c,1);           % dimension of state vector
m   = size(u0, 1);          % dimension of control vector
m_c = size(u0_c, 1);          % dimension of control vector

N   = size(u0, 2);          % number of state transitions
u   = u0;                   % initial control sequence


% --- initialize trace data structure
% trace = struct('iter',nan,'cost',nan,...
%         'dcost',nan);
% trace = repmat(trace,[min(Op.maxIter,1e6) 1]);
% trace(1).iter = 1;

iter = 1;
% --- initial trajectory for robot and contact model
[x,un,c0]         = traj_sim(x0, u0, u0_c, DYNCST);
[x_c, un_c, c0_c] = traj_sim_contact(x, x0_c,u0_c,DYNCST_CONTACT);

c01 = sum(c0(:));

u = un;

% user plotting
% Op.plotFn(x);

% constants, timers, counters
stop        = 0;
% dcost       = 0;
print_head  = 6; % print headings every print_head lines
last_head   = print_head;
verbosity   = 1;

if verbosity > 0
    fprintf('\n=========== begin ADMM ===========\n');
end

% Initialize dual variebles
rhao = [0.7,0.9];
rhao_c = [0.1,0.3];

alpha = 1.5;
% alphak = 1;
% yita = 0.999;
it2 = 30;

% Primal variables
xnew = x;
unew = u;
u_bar = ones(size(u));
x_bar = ones(size(x));
x_c_bar = ones(size(x_c));

alphak_v = ones(1,it2+1);

% Dual variables
x_lambda = zeros(size(x));
x_c_lambda = zeros(size(x_c));
u_lambda = zeros(size(u));

% ck = (1/roll)*norm(x_lambda-x_lambda2)^2 + (1/roll)*norm(u_lambda-u_lambda2)^2 + roll*norm(x_bar-x_bar2)^2 + roll*norm(u_bar-u_bar2)^2;
    
res_u = zeros(1,it2);
res_x = zeros(1,it2);
res_x_c = zeros(1,it2);
res_ulambda = zeros(1,it2);
res_xlambda = zeros(1,it2);
res_x_c_lambda = zeros(1,it2);

costcomp_robot   = zeros(1, it2);
costcomp_contact = zeros(1, it2); 

%% ADMM iteration
for i = 1:it2

    if i < 1000
    %% Original simple ADMM 
    %         %====== proximal operator to minimize to cost
    %         [xnew, unew, ~] = iLQG2(DYNCST, DYNCST2, x0, unew, roll, x_bar-x_lambda,u_bar-u_lambda, Op);
    % 
    %         %====== project operator to satisfy the constraint
    %         x_bar_old = x_bar;
    %         u_bar_old = u_bar;
    %         [x_bar, u_bar] = proj(xnew+x_lambda, unew+u_lambda, Op.lims);
    % 
    %         %====== dual variables update
    %         x_lambda = x_lambda + xnew - x_bar;
    %         u_lambda = u_lambda + unew - u_bar;
        
    %% ADMM with relaxtion 
        %====== proximal operator to minimize to cost
        % robot manipulator dynamics
        [xnew, unew, ~]   = iLQG_ADMM(DYNCST, DYNCST_primal, x0, unew, rhao, x_bar-x_lambda,u_bar-u_lambda, Op);
        
%         cnew = xcnew(1:2,:);
%         hnew = [xcnew(5,:);M*xcnew(3:4,:)];
%         hnew(1,:) = momentum(1,:)+h_lambda(1,:);
%         ufnew = [zeros(2,size(ufnew_s,2));ufnew_s([2 1],:);];
%                 
                
        % contact model, unew is common to both blocks
        [xnew_c, unew, ~] = iLQG_ADMM(DYNCST_CONTACT, DYNCYT_primal_CONTACT, x0_c, unew, rhao_c, x_c_bar-x_c_lambda, u_c_bar-u_c_lambda, Op);

        % block for cartesian constraints
        
        % Relaxtion
        xnew2 = alpha*xnew + (1-alpha)*x_bar;
        unew2 = alpha*unew + (1-alpha)*u_bar;
        
        % Relaxation of CONTACT model 
        xnew2_c = alpha*xnew_c + (1-alpha)*x_bar_c;

        %====== project operator to satisfy the constraint
        x_bar_old    = x_bar;
        x_c_bar_old  = x_c_bar;
        u_bar_old    = u_bar;
        
        [x_bar, u_bar] = proj(xnew2+x_lambda, unew2+u_lambda, Op.lims);
        [x_c_bar, ~]   = proj(xnew2_c+x_c_lambda, unew2+u_lambda, Op.lims);

        %====== dual variables update
        x_lambda = x_lambda + xnew2 - x_bar;
        u_lambda = u_lambda + unew2 - u_bar;
        
        x_c_lambda = x_c_lambda + xnew2_c - x_c_bar;
        

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
    % ====== residue ======= %
    res_u(:,i) = norm(unew - u_bar);
    res_x(:,i) = norm(xnew - x_bar);
    res_x_c(:,i) = norm(xnew_c - x_c_bar);
    
    res_ulambda(:,i) = rhao(2) * norm(u_bar - u_bar_old);
    res_xlambda(:,i) = rhao(1) * norm(x_bar - x_bar_old);
    res_x_c_lambda(:,i) = rhao_c(1) * norm(x_c_bar - x_c_bar_old);

    [~,~,cost_robot]        = traj_sim(x0, unew, DYNCST);
    [~,~,cost_contact]  = traj_sim_contact(x0_c, unew, DYNCST_CONTACT);
    
    costcomp_robot(:,i) = sum(cost_robot(:));
    costcomp_contact(:,i) = sum(cost_contact(:));
    
    % ====== varying penalty parameter ======== %
    %     if i > 5
    %         if res_u(:,i) > 10*res_ulambda(:,i)
    %             rhao(2) = 2*rhao(2);
    %             u_lambda = u_lambda/2;
    %             x_lambda = x_lambda/2;
    %         elseif res_ulambda(:,i) > 10*res_u(:,i)
    %             rhao(1) = rhao(1)/2;
    %             u_lambda = u_lambda*2;
    %             x_lambda = x_lambda*2;
    %         end
    %     end
end

figure(15)
ppp = 1:it2+1;
plot(ppp,alphak_v);

%% plot the residue
figure(10)
subplot(1,2,1)
l = 1:it2;
plot(l,res_u);
hold on;
plot(l,res_x);
plot(l,res_ulambda);
plot(l,res_xlambda);
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
plot(ii,unew);
hold on 
plot(ii,u_bar);

figure(12)
e1 = (unew(1,:)-u_bar(1,:))./unew(1,:);
e2 = (unew(2,:)-u_bar(2,:))./unew(2,:);
plot(ii,e1);
hold on 
plot(ii,e2);
hold off

[~,~,costnew]  = traj_sim(x0,unew,DYNCST);

%% ====== STEP 5: accept step (or not), expand or shrink trust region
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
cost           = costnew;
% Op.plotFn(x);
   
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


function [xnew,unew,cnew] = traj_sim(x0, u0, u_c, DYNCST)
% Generate the initial trajectory 

    n        = size(x0,1);
    m        = size(u0,1);
    N        = size(u0,2);

    xnew        = zeros(n,N);
    xnew(:,1)   = x0(:,1);
    unew        = u0;
    cnew        = zeros(1,N+1);
    
    for i = 1:N
        [xnew(:,i+1), cnew(:,i)]  = DYNCST(xnew(:,i), unew(:,i), u_c(:,i), i);
    end
    [~, cnew(:,N+1)] = DYNCST(xnew(:,N+1),nan(m,1), u_c(:,N), i);

function [xnew,unew,cnew] = traj_sim_contact(x_r, x0_c, u0_c, DYNCST_CONTACT)
% Generate the initial trajectory 

    n        = size(x0_c,1);
    m        = size(u0_c,1);
    N        = size(u0_c,2);

    xnew        = zeros(n,N);
    xnew(:,1)   = x0_c(:,1);
    unew        = u0_c;
    cnew        = zeros(1,N+1);
    for i = 1:N
        [xnew(:,i+1), cnew(:,i)]  = DYNCST_CONTACT(x_r(:,i), xnew(:,i), unew(:,i),i);
    end
    [~, cnew(:,N+1)] = DYNCST_CONTACT(x_r(:,N), xnew(:,N+1),nan(m,1),i);        
        

function [x2, u2] = proj(xnew, unew, lims)
    % Project operator(control-limit): simply clamp the control output
    N = size(unew, 2);
    m = size(unew, 1);
    n = size(xnew, 1);
    u2 = zeros(m,N);
    x2 = xnew;

    for i =1:N

        for j = 1:n
            if xnew(j,i+1) > lims(1,2)
                x2(j,i+1) = lims(1,2);
            elseif xnew(j,i+1) < lims(1,1)
                x2(j,i+1) = lims(1,1);
            else
                x2(j,i+1) = xnew(j,i+1);
            end
        end

        for k = 1:m
            if unew(k,i) > lims(2,2)
                u2(k,i) = lims(2,2);
            elseif unew(k,i) < lims(2,1)
                u2(k,i) = lims(2,1);
            else
                u2(k,i) = unew(k,i);
            end
    end
end
