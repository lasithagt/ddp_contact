
function [Slist_K, M_K] = manipulator_POE
    
        
    %% manipulator twists (kuka)
    W1 = [0 0 1]'; W2 = [0 -1 0]';
    w1 = W1; q1 = [0;0;0.2025];
    w2 = W2; q2 = [0;0;0.2025];
    w3 = W1; q3 = [0;0;0.2025];
    w4 = -W2; q4 = [0;0;0.2025+0.42];
    w5 = W1; q5 = [0;0;0.2025+0.42];
    w6 = W2; q6 = [0;0;0.2025+0.42+0.4];
    w7 = W1; q7 = [0;0;0.2025+0.42+0.4+0.126];

    g_st = [1 0 0 0;0 1 0 0;0 0 1 0.2025+0.42+0.4+0.126;0 0 0 1];
    M_K = g_st;

    h = 0;
    S1 = ScrewToAxis(q1,w1, h);
    S2 = ScrewToAxis(q2,w2, h);
    S3 = ScrewToAxis(q3,w3, h);
    S4 = ScrewToAxis(q4,w4, h);
    S5 = ScrewToAxis(q5,w5, h);
    S6 = ScrewToAxis(q6,w6, h);
    S7 = ScrewToAxis(q7,w7, h);

    S        = [S1, S2, S3, S4, S5, S6, S7];
    Slist_K  = S;
end
