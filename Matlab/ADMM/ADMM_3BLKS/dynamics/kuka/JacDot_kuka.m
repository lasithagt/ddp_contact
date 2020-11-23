function jacDot = JacDot_kuka(q, qd)
    h = 1e-10;
    q_ = qd * h + q;
    jacDot = (Jac_kuka(q_) - Jac_kuka(q))./h;
end