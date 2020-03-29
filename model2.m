% model1.M
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model2

syms DELTA ALFA BETTA G LAMBDAZ ETA MU LAMBDAP GAMA
syms c cp l lp k kp Z Zp Z1 Z1p Z2 Z2p Z3 Z3p Z4 Z4p P Pp

f1 = c + G*kp - (1-DELTA) * k - Z * k^ALFA * l^(1-ALFA);
f2 = c^(-1) - (BETTA/G) * cp^(-1) * ((1/Pp) * Zp * ALFA * kp^(ALFA-1)*lp^(1-ALFA) + 1 - DELTA);
f3 = GAMA * l^ETA * c - (1/P) * (1-ALFA) * Z * k^ALFA * l^(-ALFA);
f4 = Zp - Z^LAMBDAZ;
f5 = Pp - G * P^LAMBDAP * Z^(MU) * Z1^(MU^2) * Z2^(MU^3) * Z3^(MU^4) * Z4^(MU^5);
f6 = Z1p - Z;
f7 = Z2p - Z1;
f8 = Z3p - Z2;
f9 = Z4p - Z3;

f = [f1;f2;f3;f4;f5;f6;f7;f8;f9];

x = [k Z Z1 Z2 Z3 Z4 P];
y = [l c];
xp = [kp Zp Z1p Z2p Z3p Z4p Pp];
yp = [lp cp];

% f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);