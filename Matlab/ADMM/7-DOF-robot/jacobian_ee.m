function JAC = jacobian_ee(links, q_)
  a1 = links(1); a2 = links(2);
  q = @(x)q_(x+1);
  
  JAC(1) = -a1*sin(q(0)) - a2*sin(q(0))*cos(q(1)) - a2*sin(q(1))*cos(q(0));
  JAC(2) = -a2*sin(q(0))*cos(q(1)) - a2*sin(q(1))*cos(q(0));
  JAC(3) = a1*cos(q(0)) - a2*sin(q(0))*sin(q(1)) + a2*cos(q(0))*cos(q(1));
  JAC(4) = -a2*sin(q(0))*sin(q(1)) + a2*cos(q(0))*cos(q(1));
  JAC(5) = 0;
  JAC(6) = 0;
  JAC(7) = 0;
  JAC(8) = 0;
  JAC(9) = 0;
  JAC(10) = 0;
  JAC(11) = 1;
  JAC(12) = 1;
  JAC     = reshape(JAC,2,6)';
end