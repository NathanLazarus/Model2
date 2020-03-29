function [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,PSTAR]=model2_ss_numeric(DELTA,ALFA,BETTA,G,ETA,ZSTAR,LAMBDAP,GAMA)
% This program computes the steady state
PSTAR = G^(1/(1-LAMBDAP));
syms C K L
ss1 = C + G*K - (1-DELTA) * K - ZSTAR * K^ALFA * L^(1-ALFA);
ss2 = 1 - (BETTA/G) * ((1/PSTAR) * ZSTAR * ALFA * K^(ALFA-1)*L^(1-ALFA) + 1 - DELTA);
ss3 = GAMA * L^ETA * C - (1/PSTAR) * (1-ALFA) * ZSTAR * K^ALFA * L^(-ALFA);

[SOL]=solve(ss1,ss2,ss3,C,K,L,'Real',true)
p_alt = PSTAR;
while length(SOL.L)<1
    p_alt = p_alt+0.000001
    ss1 = C + G*K - (1-DELTA) * K - ZSTAR * K^ALFA * L^(1-ALFA);
    ss2 = 1 - (BETTA/G) * ((1/p_alt) * ZSTAR * ALFA * K^(ALFA-1)*L^(1-ALFA) + 1 - DELTA);
    ss3 = L^ETA * C - (1/p_alt) * (1-ALFA) * ZSTAR * K^ALFA * L^(-ALFA);
    [SOL]=solve(ss1,ss2,ss3,C,K,L,'Real',true)
end


CSTAR=double(SOL.C);
KSTAR=double(SOL.K);
LSTAR=double(SOL.L)
WSTAR=ZSTAR/PSTAR*(1-ALFA)*KSTAR^ALFA*LSTAR^(-ALFA)
RSTAR=ZSTAR/PSTAR*ALFA*KSTAR^(ALFA-1)*LSTAR^(1-ALFA)-DELTA;

end
