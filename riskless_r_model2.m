function riskless_r = riskless_r_model2(state_vars,c_input,l_input,dec_c,dec_l,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym)

k = state_vars(1);
Z = state_vars(2);
Z1 = state_vars(3);
Z2 = state_vars(4);
Z3 = state_vars(5);
P = state_vars(6);

quadpoints = [-0.381186990207322116,...
0.3811869902073221168,...
-1.157193712446780194,...
1.1571937124467801947,...
-1.981656756695842925,...
1.9816567566958429258,...
-2.930637420257244019,...
2.9306374202572440192];

quadweights = [0.661147012558241291,...
0.661147012558241291,...
0.207802325814891879,...
0.207802325814891879,...
0.017077983007413475,...
0.017077983007413475,...
0.000199604072211367,...
0.000199604072211367];

%points and weights from https://github.com/sivaramambikasaran/Quadrature/blob/master/Gauss_Hermite/weights/weights8
   
quad = 0;
c = c_input; %used in eval(intertemporal_euler_sym)
l = l_input; %used in the multiplicative utility case (l_plus and l don't cancel)

% for i=1:length(quadpoints)
%     Zp = Z^LAMBDAZ*exp(sqrt(2)*quadpoints(i)).^sigma_Z;
%     for j=1:length(quadpoints)
%          Pp = G * P^LAMBDAP * Zp^(MU) * Z^(MU^2) * Z1^(MU^3) * Z2^(MU^4) * Z3^(MU^5) * exp(sqrt(2)*quadpoints(j))^sigma_P;
%          cp = decision_func(dec_c,[k Zp Z Z1 Z2 Pp],ssvals,sigma_Z);
%          lp = decision_func(dec_l,[k Zp Z Z1 Z2 Pp],ssvals,sigma_Z);
%          quad = quad + quadweights(i) * (1/sqrt(pi)) * quadweights(j)*(1/sqrt(pi))*(1/(BETTA*eval(intertemporal_euler_sym)));
%     end
% %     Pp = G * P^LAMBDAP * Zp^(MU) * Z^(MU^2) * Z1^(MU^3) * Z2^(MU^4) * Z3^(MU^5);
% %     cp = decision_func(dec_c,[k Zp Z Z1 Z2 Pp],ssvals,sigma_Z);
% %     lp = decision_func(dec_l,[k Zp Z Z1 Z2 Pp],ssvals,sigma_Z);
% %    quad = quad + quadweights(i) * (1/sqrt(pi)) * (1/(BETTA*eval(intertemporal_euler_sym)));
% end
quad = 0;
for i=1:length(quadpoints)
    Zp = Z^LAMBDAZ*exp(sqrt(2)*quadpoints(i)).^sigma_Z;
    cp = decision_func(dec_c,[k Z^LAMBDAZ*exp(sqrt(2)*quadpoints(i)).^sigma_Z Z Z1 Z2 G * P^LAMBDAP * Zp^(MU) * Z^(MU^2) * Z1^(MU^3) * Z2^(MU^4) * Z3^(MU^5)],ssvals,sigma_Z);
    lp = decision_func(dec_l,[k Z^LAMBDAZ*exp(sqrt(2)*quadpoints(i)).^sigma_Z Z Z1 Z2 G * P^LAMBDAP * Zp^(MU) * Z^(MU^2) * Z1^(MU^3) * Z2^(MU^4) * Z3^(MU^5)],ssvals,sigma_Z);
    quad = quad + quadweights(i)*(1/sqrt(pi))*(1/(BETTA*eval(intertemporal_euler_sym)));
end
riskless_r = quad - 1;