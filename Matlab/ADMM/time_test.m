q_pos = ones(7,1); q_vel = 10*rand(7,1);
params = inertial_params;
tic
M    = M_kuka(params, q_pos');
CX2  = C_kuka(params, q_pos', q_vel');
Grav = -G_kuka(params, q_pos);

qdd = (M+0.000*eye(7))\( - Grav' - CX2'*q_vel)
toc

tic
M1 = massMatrix(lbr, q_pos);
C2 = velocityProduct(lbr, q_pos, q_vel);
G2 = gravityTorque(lbr, q_pos);
qdd = (M1+0.000*eye(7))\( - G2 - C2)
toc

tic
forwardDynamics(lbr, q_pos, q_vel, zeros(7,1))
toc