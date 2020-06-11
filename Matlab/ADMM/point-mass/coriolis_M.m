function C = coriolis_M(parms_, links_, q_, dq_)

  parms = @(q)parms_(q+1);
  q  = @(x)q_(x+1);
  dq = @(x)dq_(x+1);
  
  a1 = links_(1); a2 = links_(2);
  
  x0 = sin(q(1));
  x1 = a1*x0;
  x2 = cos(q(1));
  x3 = a2*x0^2 + a2*x2^2;
  x4 = a1*x2 + x3;
  x5 = parms(21)*x4;
  x6 = -parms(19);
  x7 = parms(21)*x1;
  x8 = x6 + x7;
  x9 = x1*(parms(18) + x5) - x4*x8;
  x10 = a1*parms(7);
  x11 = -parms(7);
  x12 = -parms(18);
  x13 = a1*(x0*(x12 - x5) + x11 + x2*x8) + x10 + x3*x8 + x9;
  x14 = -x1;
  x15 = -4*parms(19) + parms(21)*x14 + 2*x7;
  x16 = a2 + x4;
  x17 = parms(21)*x16;
  x18 = -a2;
  x19 = x18 + x4;
  x20 = parms(18)*x14 + x1*(2*parms(18) + x17) - x16*(-2*parms(19) + x7) + x19*x6;
  x21 = a2*parms(19);
  x22 = a1*(x0*(parms(21)*x18 + x12) + x2*x6) + x21 + x3*x6;
  x23 = a2*x8 + x9;
  x24 = a2*x6 + x21;

  C(1) = dq(0)*x13;
  C(2) = dq(0)*(a1*(x0*(-4*parms(18) + parms(21)*x19 - 2*x17) + x11 + x15*x2) + x10 - x13 + x15*x3 + x20 - x22) + dq(1)*x22;
  C(3) = dq(0)*x23; 
  C(4) = dq(0)*(a2*x15 + x20 - x23 - x24) + dq(1)*x24;
  
  C = reshape(C,2,2)';

end