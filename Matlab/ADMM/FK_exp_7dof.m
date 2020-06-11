
function rb = FK_exp_7dof
    
        
    %% manipulator twists
    W1 = [0 0 1]'; W2 = [0 -1 0];
    w1 = W1; q1 = [0;0;0.2025];
    w2 = W2; q2 = [0;0;0.2025];
    w3 = W1; q3 = [0;0;0.2025];
    w4 = -W2; q4 = [0;0;0.2025+0.42];
    w5 = W1; q5 = [0;0;0.2025+0.42];
    w6 = W2; q6 = [0;0;0.2025+0.42+0.4];
    w7 = W1; q7 = [0;0;0.2025+0.42+0.4+0.126];
    
    l1 = createtwist(w1, q1);
    l2 = createtwist(w2, q2);
    l3 = createtwist(w3, q3);
    l4 = createtwist(w4, q4);
    l5 = createtwist(w5, q5);
    l6 = createtwist(w6, q6);
    l7 = createtwist(w7, q7);

    g_st =  [1 0 0 0;0 1 0 0;0 0 1 0.2025+0.42+0.4+0.126;0 0 0 1];
    rb   =  robot({l1, l2, l3, l4, l5, l6, l7}, g_st);
    
%     q    = [0 pi/4 0 pi/10 pi/4 pi/9 0];
%     assert(fkine(rb,q')==FK_kuka(q))
%     fkine(rb,q')
%     FK_kuka(q)

end
