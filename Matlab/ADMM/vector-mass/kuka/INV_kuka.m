function INV = INV_kuka(parms_, q_, dq_, ddq_)

  q = @(n)q_(n+1);
  dq = @(n)dq_(n+1);
  ddq = @(n)ddq_(n+1);
  
  parms = @(n)parms_(n+1);
  
 x0 = sin(q(1));
  x1 = cos(q(1));
  x2 = dq(0)*x1;
  x3 = -x2;
  x4 = -dq(1);
  x5 = ddq(0)*x0 + x3*x4;
  x6 = -ddq(1);
  x7 = dq(0)*x0;
  x8 = ddq(0)*x1 + x4*x7;
  x9 = 9.81*x1;
  x10 = parms(14)*x7 + parms(16)*x4 + parms(17)*x2;
  x11 = cos(q(2));
  x12 = x11*x5;
  x13 = sin(q(2));
  x14 = x13*x6;
  x15 = x11*x4;
  x16 = x13*x7;
  x17 = x15 - x16;
  x18 = -dq(2);
  x19 = -x18;
  x20 = x12 + x14 + x17*x19;
  x21 = -ddq(2) - x8;
  x22 = x11*x6;
  x23 = x13*x5;
  x24 = x11*x7;
  x25 = x13*x4;
  x26 = x24 + x25;
  x27 = x18*x26 + x22 - x23;
  x28 = -x13;
  x29 = 9.81*x0;
  x30 = 0.42*x15 - 0.42*x16;
  x31 = -0.42*x12 - 0.42*x14 + x18*x30 + x28*x29;
  x32 = x18 + x3;
  x33 = -parms(31);
  x34 = parms(26)*x26 + parms(28)*x32 + parms(29)*x17 + x30*x33;
  x35 = cos(q(3));
  x36 = sin(q(3));
  x37 = -x35;
  x38 = x26*x36 + x32*x37;
  x39 = -dq(3);
  x40 = x20*x35 + x21*x36 + x38*x39;
  x41 = ddq(3) + x27;
  x42 = x26*x35 + x32*x36;
  x43 = dq(3)*x42 + x20*x36 + x21*x37;
  x44 = -0.42*x24 - 0.42*x25;
  x45 = x11*x29 + x19*x44 + 0.42*x22 - 0.42*x23;
  x46 = -x9;
  x47 = x30*x35;
  x48 = dq(3)*x47 + x36*x45 + x37*x46;
  x49 = x30*x36;
  x50 = -parms(42);
  x51 = dq(3) + x17;
  x52 = parms(43)*x42 + parms(45)*x49 + x50*x51;
  x53 = cos(q(4));
  x54 = x40*x53;
  x55 = sin(q(4));
  x56 = x41*x55;
  x57 = x51*x53;
  x58 = x42*x55;
  x59 = -x57 + x58;
  x60 = -dq(4);
  x61 = x54 + x56 + x59*x60;
  x62 = ddq(4) + x43;
  x63 = x40*x55;
  x64 = -x53;
  x65 = x42*x53;
  x66 = x51*x55;
  x67 = x65 + x66;
  x68 = dq(4)*x67 + x41*x64 + x63;
  x69 = x35*x45 + x36*x46 + x39*x49;
  x70 = x44*x55 + x47*x53 + 0.4*x57 - 0.4*x58;
  x71 = dq(4)*x70 + x31*x64 + 0.4*x54 + x55*x69 + 0.4*x56;
  x72 = dq(4) + x38;
  x73 = -x70;
  x74 = parms(50)*x67 + parms(52)*x72 + parms(53)*x59 + parms(54)*x49 + parms(55)*x73;
  x75 = cos(q(5));
  x76 = sin(q(5));
  x77 = -x76;
  x78 = x67*x77 + x72*x75;
  x79 = -dq(5);
  x80 = -x79;
  x81 = x61*x75 + x62*x76 + x78*x80;
  x82 = -x68;
  x83 = -ddq(5) + x82;
  x84 = x67*x75 + x72*x76;
  x85 = x61*x77 + x62*x75 + x79*x84;
  x86 = 0.4*x53;
  x87 = x44*x64 + x47*x55 + 0.4*x65 + 0.4*x66;
  x88 = x31*x55 + x41*x86 + x53*x69 + x60*x87 - 0.4*x63;
  x89 = x49*x76 + x70*x75;
  x90 = x48*x75 + x77*x88 + x79*x89;
  x91 = -x87;
  x92 = x49*x75 + x70*x77;
  x93 = -x59;
  x94 = x79 + x93;
  x95 = -parms(66);
  x96 = parms(67)*x84 + parms(69)*x92 + x94*x95;
  x97 = -parms(67);
  x98 = parms(62)*x84 + parms(64)*x94 + parms(65)*x78 + parms(66)*x91 + x89*x97;
  x99 = cos(q(6));
  x100 = x81*x99;
  x101 = sin(q(6));
  x102 = x101*x83;
  x103 = x94*x99;
  x104 = x101*x84;
  x105 = x103 - x104;
  x106 = dq(6)*x105 + x100 + x102;
  x107 = x83*x99;
  x108 = -x101;
  x109 = x84*x99;
  x110 = x101*x94;
  x111 = x109 + x110;
  x112 = -dq(6);
  x113 = x107 + x108*x81 + x111*x112;
  x114 = ddq(6) + x85;
  x115 = dq(6) + x78;
  x116 = x108*x89 - 0.126*x109 - 0.126*x110 + x91*x99;
  x117 = x101*x91 + 0.126*x103 - 0.126*x104 + x89*x99;
  x118 = -parms(79);
  x119 = parms(74)*x111 + parms(76)*x105 + parms(77)*x115 + parms(78)*x116 + x117*x118;
  x120 = -parms(78);
  x121 = parms(79)*x111 + parms(81)*x92 + x105*x120;
  x122 = -x71;
  x123 = x48*x76 + x75*x88 + x80*x92;
  x124 = -0.126*x100 - 0.126*x102 + x108*x123 + x112*x117 + x122*x99;
  x125 = -parms(80);
  x126 = parms(78)*x115 + parms(81)*x116 + x111*x125;
  x127 = -x126;
  x128 = parms(73)*x111 + parms(75)*x105 + parms(76)*x115 + parms(80)*x117 + x120*x92;
  x129 = parms(72)*x106 + parms(73)*x113 + parms(74)*x114 + parms(79)*x90 + x105*x119 - x115*x128 + x116*x121 + x124*x125 + x127*x92;
  x130 = -parms(68);
  x131 = parms(61)*x84 + parms(63)*x94 + parms(64)*x78 + parms(68)*x89 + x92*x95;
  x132 = -x78;
  x133 = parms(66)*x78 + parms(69)*x91 + x130*x84;
  x134 = -0.126*x101;
  x135 = dq(6)*x116 + x101*x122 + 0.126*x107 + x123*x99 + x134*x81;
  x136 = parms(80)*x105 + parms(81)*x117 + x115*x118;
  x137 = parms(72)*x111 + parms(73)*x105 + parms(74)*x115 + parms(79)*x92 + x116*x125;
  x138 = -x111;
  x139 = parms(73)*x106 + parms(75)*x113 + parms(76)*x114 + parms(80)*x135 + x115*x137 - x117*x121 + x119*x138 + x120*x90 + x136*x92;
  x140 = parms(80)*x113 + parms(81)*x135 + x105*x121 + x114*x118 + x115*x127;
  x141 = parms(78)*x114 + parms(81)*x124 + x106*x125 + x115*x136 + x121*x138;
  x142 = x141*x99;
  x143 = parms(60)*x81 + parms(61)*x83 + parms(62)*x85 + parms(67)*x90 + x108*x139 + x122*x130 + x129*x99 + x131*x132 - x133*x92 + x134*x140 - 0.126*x142 + x91*x96 + x94*x98;
  x144 = -x72;
  x145 = parms(54)*x144 + parms(55)*x67 + parms(57)*x87;
  x146 = -parms(56);
  x147 = parms(49)*x67 + parms(51)*x72 + parms(52)*x59 + parms(54)*x91 + parms(56)*x70;
  x148 = parms(54)*x59 + parms(57)*x49 + x146*x67;
  x149 = -x105;
  x150 = parms(74)*x106 + parms(76)*x113 + parms(77)*x114 + parms(78)*x124 + x111*x128 - x116*x136 + x117*x126 + x118*x135 + x137*x149;
  x151 = parms(68)*x94 + parms(69)*x89 + x78*x97;
  x152 = -x151;
  x153 = parms(60)*x84 + parms(61)*x94 + parms(62)*x78 + parms(67)*x92 + x130*x91;
  x154 = parms(62)*x81 + parms(64)*x83 + parms(65)*x85 + parms(66)*x122 + x123*x97 + x131*x84 + x133*x89 + x150 + x152*x91 - x153*x94;
  x155 = parms(48)*x61 + parms(49)*x62 + parms(50)*x68 + parms(55)*x71 + x143*x75 + x145*x49 + x146*x48 + x147*x93 + x148*x91 + x154*x77 + x72*x74;
  x156 = -x84;
  x157 = x140*x99;
  x158 = -parms(61)*x81 - parms(63)*x83 - parms(64)*x85 - parms(68)*x123 - x101*x129 - x134*x141 - x139*x99 - x151*x92 - x153*x78 - x156*x98 - 0.126*x157 + x89*x96 - x90*x95;
  x159 = parms(48)*x67 + parms(49)*x72 + parms(50)*x59 + parms(55)*x87 + x146*x49;
  x160 = parms(55)*x93 + parms(56)*x72 + parms(57)*x70;
  x161 = -x49;
  x162 = parms(50)*x61 + parms(52)*x62 + parms(53)*x68 + parms(54)*x48 - parms(55)*x88 + x144*x159 + x147*x67 + x148*x70 + x158 + x160*x161;
  x163 = -parms(43);
  x164 = parms(38)*x42 + parms(40)*x51 + parms(41)*x38 + parms(42)*x44 + x163*x47;
  x165 = -parms(44);
  x166 = parms(37)*x42 + parms(39)*x51 + parms(40)*x38 + parms(42)*x161 + parms(44)*x47;
  x167 = -x38;
  x168 = parms(42)*x38 + parms(45)*x44 + x165*x42;
  x169 = -parms(54)*x62 + parms(55)*x61 + parms(57)*x71 - parms(66)*x85 - parms(69)*x122 - x101*x140 - x130*x81 - x142 + x144*x160 + x148*x67 - x151*x78 - x156*x96;
  x170 = parms(68)*x83 + parms(69)*x123 + x108*x141 + x132*x133 + x157 + x85*x97 + x94*x96;
  x171 = parms(67)*x81 + parms(69)*x90 + parms(79)*x106 + parms(81)*x90 + x111*x126 + x113*x120 + x133*x84 + x136*x149 + x152*x94 + x83*x95;
  x172 = parms(55)*x82 + parms(56)*x62 + parms(57)*x88 + x145*x72 + x148*x93 + x170*x75 + x171*x77;
  x173 = x172*x55;
  x174 = parms(36)*x40 + parms(37)*x41 + parms(38)*x43 + parms(43)*x48 + x155*x53 + x161*x168 + x162*x55 + x164*x51 + x165*x31 + x166*x167 + x169*x86 - 0.4*x173 + x44*x52;
  x175 = -x67;
  x176 = parms(49)*x61 + parms(51)*x62 + parms(52)*x68 + parms(54)*x122 + parms(56)*x88 + x143*x76 + x145*x73 + x154*x75 + x159*x59 + x160*x87 + x175*x74;
  x177 = parms(44)*x51 + parms(45)*x47 + x163*x38;
  x178 = -x177;
  x179 = parms(36)*x42 + parms(37)*x51 + parms(38)*x38 + parms(43)*x49 + x165*x44;
  x180 = parms(38)*x40 + parms(40)*x41 + parms(41)*x43 + parms(42)*x31 + x163*x69 + x166*x42 + x168*x47 + x176 + x178*x44 - x179*x51;
  x181 = -parms(32);
  x182 = -x44;
  x183 = parms(25)*x26 + parms(27)*x32 + parms(28)*x17 + parms(30)*x182 + parms(32)*x30;
  x184 = -x17;
  x185 = parms(30)*x17 + x181*x26;
  x186 = parms(24)*x20 + parms(25)*x21 + parms(26)*x27 + parms(31)*x31 + x174*x35 + x180*x36 + x181*x46 + x182*x185 + x183*x184 + x32*x34;
  x187 = parms(13)*x7 + parms(15)*x4 + parms(16)*x2;
  x188 = -x42;
  x189 = x169*x55;
  x190 = parms(37)*x40 + parms(39)*x41 + parms(40)*x43 + parms(44)*x69 + x155*x55 + x162*x64 + x164*x188 + x172*x86 + x177*x49 + x179*x38 + 0.4*x189 - x47*x52 + x48*x50;
  x191 = parms(24)*x26 + parms(25)*x32 + parms(26)*x17 + parms(31)*x44;
  x192 = -x32;
  x193 = parms(26)*x20 + parms(28)*x21 + parms(29)*x27 + parms(30)*x46 + x183*x26 + x185*x30 + x190 + x191*x192 + x33*x45;
  x194 = -parms(30);
  x195 = parms(31)*x26 + parms(33)*x44 + x194*x32;
  x196 = parms(32)*x21 + parms(33)*x45 + x184*x185 + x195*x32 + x27*x33 + x35*(parms(44)*x41 + parms(45)*x69 + x163*x43 + x167*x168 + x172*x53 + x189 + x51*x52) + x36*(parms(43)*x40 + parms(45)*x48 + parms(54)*x68 + parms(57)*x48 + x145*x175 + x146*x61 + x160*x59 + x168*x42 + x170*x76 + x171*x75 + x178*x51 + x41*x50);
  x197 = -0.42*x13;
  x198 = parms(32)*x32 + parms(33)*x30 + x17*x33;
  x199 = parms(31)*x20 + parms(33)*x31 + parms(42)*x43 + parms(45)*x31 + x165*x40 + x169*x64 + x173 + x177*x38 + x185*x26 + x188*x52 + x192*x198 + x194*x21;
  x200 = -parms(25)*x20 - parms(27)*x21 - parms(28)*x27 - parms(32)*x45 - x17*x191 - x174*x36 - x180*x37 - x194*x31 + x195*x30 - x198*x44 + x26*x34;
  x201 = parms(12)*x7 + parms(13)*x4 + parms(14)*x2;

  INV(1) = ddq(0)*parms(3) + dq(0)*parms(10) + parms(11)*sign(dq(0)) + x0*(parms(12)*x5 + parms(13)*x6 + parms(14)*x8 + parms(19)*x9 + x10*x4 + x11*x186 - 0.42*x11*x199 + x187*x3 + x193*x28 + x196*x197) + x1*(parms(14)*x5 + parms(16)*x6 + parms(17)*x8 - parms(19)*x29 + x187*x7 + x200 - x201*x4);
  INV(2) = dq(1)*parms(22) - parms(13)*x5 - parms(15)*x6 - parms(16)*x8 - parms(18)*x46 - parms(20)*x29 + parms(23)*sign(dq(1)) + x10*x7 - x11*x193 - 0.42*x11*x196 - x13*x186 - x197*x199 - x2*x201;
  INV(3) = dq(2)*parms(34) + parms(35)*sign(dq(2)) + x200;
  INV(4) = dq(3)*parms(46) + parms(47)*sign(dq(3)) + x190;
  INV(5) = dq(4)*parms(58) + parms(59)*sign(dq(4)) + x176;
  INV(6) = dq(5)*parms(70) + parms(71)*sign(dq(5)) + x158;
  INV(7) = dq(6)*parms(82) + parms(83)*sign(dq(6)) + x150;


end