T    = 400;
t    = linspace(0, 2*pi, T+1);
r    = 0.06;
xd_x = r * cos(t);
xd_y = r * sin(2*t);
xd_z = (0.8) * ones(1,numel(t));
Rd_r = 0 * ones(1, numel(t));
Rd_p = 0 * ones(1, numel(t));
Rd_y = 0 * ones(1, numel(t));

x_des = [xd_x; xd_y; xd_z; Rd_r; Rd_p; Rd_y];

T0          = [[eye(3);0 0 0],[xd_x(1) xd_y(1) xd_z(1) 1]'];
theta0      = 0 + ones(7,1)*0.1;
theta0(2)   = 0.2;
theta0(4)   = 0.2;
[Slist, M_] = manipulator_POE();
[q0,s]      = IKinSpace_modified_initial(Slist, M_, T0, theta0, 0.00001, 0.00001);

q_bar  = zeros(7, size(x_des,2));
qd_bar = zeros(7, size(x_des,2));

bars = load('bar_data.mat');
q_bar = bars.q_bar;
qd_bar = bars.qd_bar;

[q_bar, qd_bar, ~]  = kuka_second_order_IK(x_des, q0, zeros(7,1), [0;1.6], q_bar, qd_bar, 1);
% save('bar_data.mat','q_bar','qd_bar')