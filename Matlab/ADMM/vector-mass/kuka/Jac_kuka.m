function jac = Jac_kuka(q_)

  q = @(n)q_(n+1);

  jac(1) = 0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) - 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) - 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) + 0.4*sin(q(0))*sin(q(1))*cos(q(3)) + 0.42*sin(q(0))*sin(q(1));
  jac(2) = -(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3)) + 0.42*cos(q(1)))*cos(q(0));
  jac(3) = -(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)) - 0.42*sin(q(0))*sin(q(1)))*cos(q(1)) - (-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3)) + 0.42*cos(q(1)))*sin(q(0))*sin(q(1));
  jac(4) = (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3))) + (-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)))*sin(q(1))*sin(q(2));
  jac(5) = ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3))) - (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)));
  jac(6) = -(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)))*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*sin(q(4)) + sin(q(1))*sin(q(2))*cos(q(4))) + (((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*sin(q(4)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*cos(q(4)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)));
  jac(7) = (-(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5))) - (-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)))*(-((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)));
  jac(8) = -0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)) - 0.42*sin(q(1))*cos(q(0));
  jac(9) = -(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3)) + 0.42*cos(q(1)))*sin(q(0));
  jac(10) = (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)) - 0.42*sin(q(1))*cos(q(0)))*cos(q(1)) + (-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3)) + 0.42*cos(q(1)))*sin(q(1))*cos(q(0));
  jac(11) = -(-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3))) - (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)))*sin(q(1))*sin(q(2));
  jac(12) = -((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3))) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*(-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)));
  jac(13) = (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*sin(q(4)) + sin(q(1))*sin(q(2))*cos(q(4))) - (((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*sin(q(4)) - (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(4)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)));
  jac(14) = -(-(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*(-0.257*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.257*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5))) + (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*(-((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)));
  jac(15) = 0;
  jac(16) = (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)) - 0.42*sin(q(1))*cos(q(0)))*cos(q(0)) + (-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)) - 0.42*sin(q(0))*sin(q(1)))*sin(q(0));
  jac(17) = (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)) - 0.42*sin(q(1))*cos(q(0)))*sin(q(0))*sin(q(1)) - (-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)) - 0.42*sin(q(0))*sin(q(1)))*sin(q(1))*cos(q(0));
  jac(18) = (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3))) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*(-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)));
  jac(19) = ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3))) - ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*(-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)));
  jac(20) = -(-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*sin(q(4)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*cos(q(4))) + (-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)))*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*sin(q(4)) - (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(4)));
  jac(21) = (-(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*(-0.257*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.257*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5))) - (-0.257*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.257*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)))*(-(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)));
  jac(22) = 0;
  jac(23) = sin(q(0));
  jac(24) = -sin(q(1))*cos(q(0));
  jac(25) = -sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1));
  jac(26) = (-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3));
  jac(27) = ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*sin(q(4)) - (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(4));
  jac(28) = -(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5));
  jac(29) = 0;
  jac(30) = -cos(q(0));
  jac(31) = -sin(q(0))*sin(q(1));
  jac(32) = -sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2));
  jac(33) = (sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3));
  jac(34) = ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*sin(q(4)) - (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*cos(q(4));
  jac(35) = -(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5));
  jac(36) = 1;
  jac(37) = 0;
  jac(38) = cos(q(1));
  jac(39) = -sin(q(1))*sin(q(2));
  jac(40) = sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3));
  jac(41) = (sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*sin(q(4)) + sin(q(1))*sin(q(2))*cos(q(4));
  jac(42) = -((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5));

  jac = reshape(jac,6,[]);
end