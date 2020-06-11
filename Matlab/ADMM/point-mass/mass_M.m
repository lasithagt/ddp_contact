function M = mass_M(parms_, links, q_)
    
  parms = @(x)parms_(x+1);
  q = @(x)q_(x+1);
  
  a1 = links(1);
  a2 = links(2);
  
  x0 = cos(q(1));
  x1 = sin(q(1));
  x2 = a2*x0^2 + a2*x1^2;
  x3 = a1*x0 + x2;
  x4 = a1*x1;
  x5 = -parms(19);
  x6 = parms(17) + parms(18)*x3 + x4*x5;
  x7 = parms(18) + parms(21)*x3;
  x8 = a2*x7 + x6;

  M(1) = a1*parms(6) + a1*(a1*parms(9) + parms(6) + x0*x7 + x1*(parms(21)*x4 + x5)) + parms(5) + x2*x7 + x6;
  M(2) = x8;
  M(3) = x8;
  M(4) = a2*parms(18) + a2*(a2*parms(21) + parms(18)) + parms(17);

  M = reshape(M,2,2);
end