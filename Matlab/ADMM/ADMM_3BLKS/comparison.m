clear;
T       = 400;              % horizon 

cu_w  = 5e-10*ones(1,7);  % control cost coefficients
cx_w  = 5e-1 * [0.0*ones(1,7) 0.0001*ones(1,7) 0.000000 0.00000 0.5];     % running cost coefficients


% desired path to track
t    = linspace(0, 2*pi, T+1);
r    = 0.06;
xd_x = r * cos(t);
xd_y = r * sin(3*t);
xd_z = (0.8) * ones(1,numel(t));
Rd_r = 0 * ones(1, numel(t));
Rd_p = 0 * ones(1, numel(t));
Rd_y = 0 * ones(1, numel(t));
% [xd_x, xd_y, xd_z] = lissajous_curve(t, 1.1608);

xd_f = -5.0 * ones(1,numel(t));

x_des = [xd_x; xd_y; xd_z; Rd_r; Rd_p; Rd_y;xd_f];

Op.x_des = x_des;
Op.T     = T;
Op.lims  = [0 2*pi;     % position
            -0.1 0.1;   % velocity
            -10 10;     % force
             -3  3];    % acceleration limits (m/s^2)
Op.plot  = 1;           % plot the derivatives as well
Op.maxIter = 10;
Op.cx_w = cx_w;
Op.cu_w = cu_w;

mode = 1;
switch mode

    case 1
        % admm three blocks consensus
        Op.rhao   = [3, 1e-1, 0, 0, 3];
        Op.admmMaxIter  = 10;

        demo_dynamics_ADMM_3BLKS(Op)
        fprintf('\nADMM 3 BLOCKS CONSENSUS\n')
    case 2

        % admm three blocks sequential
        Op.rhao   = [1e-3, 1, 2, 0, 2];
        Op.admmMaxIter  = 10;

        demo_dynamics_ADMM_3BLKS_SEQUENTIAL(Op)
        fprintf('\nADMM 3 BLOCKS SEQUENTIAL\n')
    case 3
        % admm two blocks 
        Op.cu_w   = [5e-10 * ones(1, 7)];
        Op.cx_w   = 5e-1 * [0.0*ones(1,7) 0.0001*ones(1,7) 0.000000 0.0001 5];
        Op.rhao   = [1e-3, 1e-1, 1, 0, 4];
        Op.admmMaxIter  = 20;

        demo_dynamics_ADMM_2BLKS(Op)
        fprintf('\nADMM 2 BLOCKS SEQENNTIAL\n')
    case 4
        % admm two blocks sequential dynamic time
        Op.cu_w   = [5e-10 * ones(1, 7) 0.4];
        Op.cx_w   = 5e-1 * [0.0*ones(1,7) 0.0001*ones(1,7) 0.000000 0.0001 5 0.001];
        Op.rhao   = [1e-3, 1e-1, 3, 0, 4];
        Op.admmMaxIter  = 20;

        demo_dynamics_ADMM_2BLKS_DYN_TIME(Op)
        fprintf('\nADMM 2 BLOCKS SEQENNTIAL DYN TIME\n')
    case 5
        % admm three blocks consensus dynamic time
        Op.cu_w   = [5e-10 * ones(1, 7) 0.4];
        Op.cx_w   = 5e-1 * [0.0*ones(1,7) 0.0001*ones(1,7) 0.000000 0.0001 5 0.001];
        Op.rhao   = [1e-3, 1e-1, 3, 0, 2];
        Op.admmMaxIter  = 20;

        demo_dynamics_ADMM_3BLKS_DYN(Op)
        fprintf('\nADMM 3 BLOCKS CONSENSUS DYN TIME\n')
end
