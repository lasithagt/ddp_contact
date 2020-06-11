close all
figure(1)

tool = [0 0 0.136]';
[Slist, M] = manipulator_POE();
theta0 = 0 + abs(0.1*rand(7,1));
theta0_ = theta0;
theta0(2) = 0.05;
theta0(3) = 0.1;
theta0(4) = 0.1;
theta0(6) = 0.1;
fileID_q  = fopen('./desired_position_wrench_.txt','w');

T_b = [[rotx(pi);[0 0 0]],[0 0 0.904+0.136 1]'];
for i = 1:size(x,2)-1
    
    R      = eul2rotm(x(4:6,i)', 'ZYX');
    v      = R' * tool;
    
    %     x_(3,i) = x(3,i) + 0.75;
    p      = x(1:3,i) + v;
    
    T          = T_b\[[R;[0 0 0]],[p;1]];
    T(1:3,1:3) = R;
    
    [q(:,i),s] = IKinSpace_modified(Slist, M, T,theta0, 0.001, 0.001);
    i
    s
    q(:,i) = mod(q(:,i),2*pi);
    
    
    mask = q(:,i)>(pi-0.1);
    q(mask,i) =  q(mask,i) - 2*pi;
%     q(mask,i) =  pi-0.1 ;
    theta0 = q(:,i);
    
    fprintf(fileID_q, '%3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n', [q(:,i);u(:,i)]);

end
