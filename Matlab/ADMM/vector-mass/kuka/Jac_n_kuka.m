function jac_w_n = Jac_n_kuka(q_)

    q = @(n)q_(n+1);

    jac_w_n(1) = -(-(-(-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(0)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(0)))*cos(q(1)) - (sin(q(0))^2*sin(q(1)) + sin(q(1))*cos(q(0))^2)*sin(q(1))*sin(q(2)))/(sin(q(0))^2*sin(q(1)) + sin(q(1))*cos(q(0))^2);
    jac_w_n(2) = ((-(-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(0)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(0)))*sin(q(0))*sin(q(1)) + (sin(q(0))^2*sin(q(1)) + sin(q(1))*cos(q(0))^2)*(-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2))))/((sin(q(0))^2*sin(q(1)) + sin(q(1))*cos(q(0))^2)*cos(q(0)));
    jac_w_n(3) = -(-(-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(0)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(0)))/(sin(q(0))^2*sin(q(1)) + sin(q(1))*cos(q(0))^2);
    jac_w_n(4) = 1;
    jac_w_n(5) = 0;
    jac_w_n(6) = 0;
    jac_w_n(7) = 0;

end