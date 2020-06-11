function FK = FK_kuka(q_)
tic
   q = @(n)q_(n+1);


  FK(1) = ((((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*cos(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*sin(q(5)))*cos(q(6)) + (-((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*sin(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(4)))*sin(q(6));
  FK(2) = -((((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*cos(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*sin(q(5)))*sin(q(6)) + (-((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*sin(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*cos(q(4)))*cos(q(6));
  FK(3) = -(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + ((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5));
  FK(4) = -0.126*(((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*cos(q(3)) + sin(q(1))*sin(q(3))*cos(q(0)))*cos(q(4)) + (-sin(q(0))*cos(q(2)) - sin(q(2))*cos(q(0))*cos(q(1)))*sin(q(4)))*sin(q(5)) + 0.126*((-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - sin(q(1))*cos(q(0))*cos(q(3)))*cos(q(5)) + 0.4*(-sin(q(0))*sin(q(2)) + cos(q(0))*cos(q(1))*cos(q(2)))*sin(q(3)) - 0.4*sin(q(1))*cos(q(0))*cos(q(3)) - 0.42*sin(q(1))*cos(q(0));
  FK(5) = ((((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*cos(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*sin(q(5)))*cos(q(6)) + (-((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*sin(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*cos(q(4)))*sin(q(6));
  FK(6) = -((((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*cos(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*sin(q(5)))*sin(q(6)) + (-((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*sin(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*cos(q(4)))*cos(q(6));
  FK(7) = -(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + ((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5));
  FK(8) = -0.126*(((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*cos(q(3)) + sin(q(0))*sin(q(1))*sin(q(3)))*cos(q(4)) + (-sin(q(0))*sin(q(2))*cos(q(1)) + cos(q(0))*cos(q(2)))*sin(q(4)))*sin(q(5)) + 0.126*((sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - sin(q(0))*sin(q(1))*cos(q(3)))*cos(q(5)) + 0.4*(sin(q(0))*cos(q(1))*cos(q(2)) + sin(q(2))*cos(q(0)))*sin(q(3)) - 0.4*sin(q(0))*sin(q(1))*cos(q(3)) - 0.42*sin(q(0))*sin(q(1));
  FK(9) = (((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*cos(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*sin(q(5)))*cos(q(6)) + (-(sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*sin(q(4)) - sin(q(1))*sin(q(2))*cos(q(4)))*sin(q(6));
  FK(10) = -(((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*cos(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*sin(q(5)))*sin(q(6)) + (-(sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*sin(q(4)) - sin(q(1))*sin(q(2))*cos(q(4)))*cos(q(6));
  FK(11) = -((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + (sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5));
  FK(12) = -0.126*((sin(q(1))*cos(q(2))*cos(q(3)) - sin(q(3))*cos(q(1)))*cos(q(4)) - sin(q(1))*sin(q(2))*sin(q(4)))*sin(q(5)) + 0.126*(sin(q(1))*sin(q(3))*cos(q(2)) + cos(q(1))*cos(q(3)))*cos(q(5)) + 0.4*sin(q(1))*sin(q(3))*cos(q(2)) + 0.4*cos(q(1))*cos(q(3)) + 0.42*cos(q(1)) + 0.2025;
  FK(13) = 0;
  FK(14) = 0;
  FK(15) = 0;
  FK(16) = 1;

   
  FK = reshape(FK,4,[])';
toc
end