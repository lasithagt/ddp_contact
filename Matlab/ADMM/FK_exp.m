
function rb = FK_exp
    
        
    %% manipulator twists
    W = [0 0 1]';
    w1 = W; q1 = [0;0;0];
    w2 = W; q2 = [0.5;0;0];
    
    %     w3 = R * [0;0;1];  q3 = Q + R * [0;0;d1];
    %     w4 = R * [0;1;0];  q4 = Q + R * [0;0;d2];
    %     w5 = R * [0;1;0];  q5 = Q + R * [0;0;d3];
    %     w6 = R * [-1;0;0]; q6 = Q + R * [0;0;d3];
    %     w7 = R * [0;0;1]; q7 = Q + R * [0;0;d3+input.tool];
    
    l1 = createtwist(w1, q1);
    l2 = createtwist(w2, q2);
    
    %     l3 = createtwist(w3, q3);
    %     l4 = createtwist(w4, q4);
    %     l5 = createtwist(w5, q5);
    %     l6 = createtwist(w6, q6);
    %     l7 = createtwist(w7, q7);

    M1 =  [1 0 0 1;0 1 0 0;0 0 1 0;0 0 0 1];
    rb =  robot({l1, l2}, M1);
    
    %     hn = fkine(rn,q)
    %     guess = zeros(7,1) + 0.1*rand(7,1);
    %     pose = eye(4);
    %     pose(3,4) = 10;
    %     theta_hat  = ikine(rn, pose, guess)
    % %     theta_hat = ikine2(rn, hn, guess)
    %     hn_ = fkine(rn,theta_hat);
    % %   
    %     x = [0,0,0.5,0.5];
    %     FK_test(q,x)
    %     if (isequalf(FK_test(q,x),hn,1e-3))
    %         fprintf(' correct solution found\n')
    %     end
end
