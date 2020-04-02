function P = P_func(G,P1,LAMBDAP,Z1,Z2,Z3,Z4,Z5,MU,rho_P)
P = G * P1^LAMBDAP * Z1^(MU) * Z2^(MU^2) * Z3^(MU^3) * Z4^(MU^4) * Z5^(MU^5)*exp(rho_P);