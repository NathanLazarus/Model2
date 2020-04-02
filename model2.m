% model1.M
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u)

syms DELTA ALFA BETTA G LAMBDAZ ETA MU LAMBDAP GAMA SIGM
syms c cp l lp k kp Z Zp Z1 Z1p Z2 Z2p Z3 Z3p Z4 Z4p P Pp

if ~exist('u','var')
    u = (1/(1 - SIGM)) * c^(1-SIGM) - GAMA/(1 + ETA) * l^(1 + ETA);
end

up = subs(u,[c l],[cp lp]);
dupdcp = jacobian(up,cp);
dudc = jacobian(u,c);
%dudl = jacobian(u,l);

f1 = c + G*kp - (1-DELTA) * k - y_func(k,l,Z,ALFA);
f2 = simplify(dudc) - (BETTA/G) * simplify(dupdcp) * big_R(kp,lp,Pp,Zp,ALFA,DELTA);
f3 = laborsupply(u) + w_func(k,l,P,Z,ALFA);
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

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);