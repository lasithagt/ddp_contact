close all
figure(1)
T = eye(4);

A = [0 0 0];
B = [0.01 0 0];
C = [0 0.01 0];
D = [0 0 0.2];

E = [0 0.01 0.2];
F = [0.01 0 0.2];
G = [0.01 0.01 0];
H = [0.01 0.01 0.2];
P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
% plot3(P(:,1),P(:,2),P(:,3),'k'), hold on % original cube
axis([-0.25 0.25 -0.25 0.25 -0.25 0.25])

roll = -0.3064; pitch = -1.2258; yaw = 9.8066;
dcm = angle2dcm(yaw, pitch, roll);
% P = P*dcm;
% plot3(P(:,1),P(:,2),P(:,3),'k') % rotated cube
axis([-0.25 0.25 -0.25 0.25 -0.25 0.25])
tool = [0 0 0.136]'/2;

% patch([1 -1 -1 1], [1 1 -1 -1], [0 0 0 0], [0.9 0.9 0.9 0.9])


hold on
for i = 1:size(x,2)
    if (mod(i,10) == 0)
        
        R = eul2rotm(x(4:6,i)', 'ZYX');
        v = R' * tool;
        % p = x(1:3,i) + v;
        p = x(7:9,i);
        hold on
        plot3(p(1),p(2),p(3),'k*')
        
        axang = [cross([0 0 1]',x(19:21,i))' acos(dot([0 0 1]',x(19:21,i))/norm(x(19:21,i))) ];
        T(1:3,1:3) = axang2rotm(axang);
        T(1:3,end) = x(1:3,i);
        % quiver3(x(1,i), x(2,i), x(3,i), x(19,i), x(20,i), x(21,i), 0.1, 'LineWidth', 2, 'Color','k');
        % hold on
        % drawframetraj(T,0.04)
        plot3(x(1,i), x(2,i), x(3,i), '*')
        axis([-0.25 0.25 -0.25 0.25 -0.25 0.25])
        
        P1 = P*R + [x(1:3,i)]';
        plot3(P1(:,1),P1(:,2),P1(:,3),'k') % rotated cube
        
        pause(0.1)
        hold on
        
    end
end
