function INV = inv_dyn_T(parms, q, dq, ddq )

  x0 = ddq(0) + ddq(1);
  x1 = cos(q(1));
  x2 = sin(q(1));
  x3 = a2*x1^2 + a2*x2^2;
  x4 = a1*ddq(0);
  x5 = a1*dq(0);
  x6 = x2*x5;
  x7 = a2*ddq(1) + ddq(0)*x3 - dq(1)*x6 + x1*x4;
  x8 = dq(0) + dq(1);
  x9 = -x8;
  x10 = parms(19)*x9 + parms(21)*x6;
  x11 = parms(18)*x0 + parms(21)*x7 + x10*x8;
  x12 = dq(0)*x3 + x1*x5;
  x13 = a2*dq(1);
  x14 = x12 + x13;
  x15 = parms(18)*x8 + parms(21)*x14;
  x16 = -dq(0)*x13 + dq(1)*x12 + x2*x4;
  x17 = -parms(19);
  x18 = parms(17)*x0 + parms(18)*x7 - x10*x14 + x15*x6 + x16*x17;

  INV(1) = a1*(ddq(0)*parms(6) - dq(0)^2*parms(7) + parms(9)*x4 + x1*x11 + x2*(parms(21)*x16 + x0*x17 + x15*x9)) + ddq(0)*parms(5) + dq(0)*parms(10) + dq(0)*parms(7)*x5 + parms(11)*sign(dq(0)) + parms(6)*x4 + x11*x3 + x18;
  INV(2) = a2*x11 + dq(1)*parms(22) + parms(23)*sign(dq(1)) + x18;

end