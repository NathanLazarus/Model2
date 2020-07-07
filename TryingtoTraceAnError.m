% model2_run.M
% Calls: model2.m num_eval.m  model2_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
% function [rgwkcl_mat] = model2_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,randomseq,order)
% 
defaults = [0.08,0.32,0.98,1.014,10,0.95,0.909,0.0109,0.02,0.37,0.5,0.3,10000,0.0109,1,0,0,2,2];
var={'DELTA','ALFA','BETTA','G','SIGM','LAMBDAP','LAMBDAZ','sigma_Z','sigma_P','MU','FRISCHELAS','STEADYSTATEL','T','shock','k0_mult','MultiplicativeU','startopposite','randomseq','order'};

DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
LAMBDAP = 0.95;
LAMBDAZ = 0.909;
sigma_Z = 0.0109;
sigma_P = 0.02;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU      = 0.37;
T       = 1000;
shock   = 0.0109;
k0_mult = 1;
MultiplicativeU = 0;
startopposite = 0;
randomseq = 2;
order = 2;

for i = 1:length(defaults)
    if ~exist(var{i},'var')
        eval(sprintf('%s = %g;',var{i},defaults(i)))
    end
end

if MultiplicativeU
    u = multiplicative_u;
else
    u = additive_u;
end

order = 2;
decision_func_to_use = @decision_func_SchmittGrohe;

LAMBDAPhigh = 0.95;
LAMBDAPlow = 0.8;
if LAMBDAP == LAMBDAPhigh
    LAMBDAPopposite = LAMBDAPlow;
end
if LAMBDAP == LAMBDAPlow
    LAMBDAPopposite = LAMBDAPhigh;
end
regimechanges = 0;
startopposite = 0;

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model2(u);

eta     = [0 0; 1 0; 0 0; 0 0; 0 0; 0 sigma_P/sigma_Z]; %Matrix defining driving force

ZSTAR = 1; %steady-state value of technology shock 
PSTAR = G^(1/(1-LAMBDAP)); %steady state markup

sym_labor_supply = laborsupply(u);
intertemporal_euler_ss = dupdcp_over_dudc(u,1);
intertemporal_euler_sym = dupdcp_over_dudc(u,0);

% addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
% import casadi.*

value_of_P_where_LSTAR_equals_STEADYSTATEL = G^(1/(1-LAMBDAPlow));

% [~,~,~,~,~,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u);
GAMA = 28.1677381532213;
% [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model2_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTAR,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss)
KSTAR = 0.835777225015446;
LSTAR = 0.273094921033578;
CSTAR = 0.312063920996547;
RSTAR = little_r(KSTAR,LSTAR,PSTAR,1,ALFA,DELTA);
% [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,PSTAR,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u)
if startopposite
    Popposite = G^(1/(1-LAMBDAPopposite));
    KSTARopposite = model2_ss_numeric(KSTAR,CSTAR,LSTAR,DELTA,ALFA,BETTA,G,Popposite,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
    k0_mult = KSTARopposite/KSTAR;
end

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR; Z1 = ZSTAR; Z2 = ZSTAR; Z3 = ZSTAR; P = PSTAR;
kp=k; cp=c; lp=l; Zp=Z; Z1p = Z; Z2p = Z; Z3p = Z; Pp = P; riskless_r_ = RSTAR;
stockSTAR = ((PSTAR-1)/PSTAR)*y_func(KSTAR,LSTAR,ZSTAR,ALFA)/(1 - BETTA*G^(1-SIGM)); stock = stockSTAR; stockp = stockSTAR;
risky_r_ = RSTAR; risky_rp_ = RSTAR;

%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

flatten = @(A) A(:);
[nstate,~] = size(hx);

if approx == 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
    
    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
    
    dec_k=[KSTAR,hx(1,:),1/2*flatten(hxx(1,:,:))',1/2*hss(1)];
    dec_l=[LSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
    dec_c=[CSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
    dec_riskless_r=[RSTAR,gx(3,:),1/2*flatten(gxx(3,:,:))',1/2*gss(3)];
    dec_stock = [stockSTAR,gx(4,:),1/2*flatten(gxx(4,:,:))',1/2*gss(4)];
    dec_risky_r = [RSTAR,gx(4,:),1/2*flatten(gxx(4,:,:))',1/2*gss(4)];

else
%     dec_k=[KSTAR,hx(1,:),zeros([1 nstate^2+1])]; 
%     dec_l=[LSTAR,gx(1,:),zeros([1 nstate^2+1])];
%     dec_c=[CSTAR,gx(2,:),zeros([1 nstate^2+1])];
    
end

% firstorder = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAP*100)),"coefs1.csv"],""));
% secondorder = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAP*100)),"coefs2.csv"],""));
% thirdorder = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAP*100)),"coefs3.csv"],""));

% if order >= 4
%     fourthorder = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAP*100)),"coefs4.csv"],""));
% end
% 
% if order == 1
%     secondorder = zeros([3 64]);
% end
% if order <= 2
%     thirdorder = zeros([3 512]);
% end
% if order <= 3
%     fourthorder = zeros([3 4096]);
% end

% dec_k=[KSTAR,firstorder(1,:),1/2*secondorder(1,:),1/6*thirdorder(1,:),1/24*fourthorder(1,:)];
% dec_l=[LSTAR,firstorder(2,:),1/2*secondorder(2,:),1/6*thirdorder(2,:),1/24*fourthorder(2,:)];
% dec_c=[CSTAR,firstorder(3,:),1/2*secondorder(3,:),1/6*thirdorder(3,:),1/24*fourthorder(3,:)];


% firstorder_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAPopposite*100)),"coefs1.csv"],""));
% secondorder_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAPopposite*100)),"coefs2.csv"],""));
% thirdorder_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAPopposite*100)),"coefs3.csv"],""));
% fourthorder_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Mu",string(round(MU*100)),"LambdaP",string(round(LAMBDAPopposite*100)),"coefs4.csv"],""));
% dec_k_opposite=[KSTAR,firstorder_opposite(1,:),1/2*secondorder_opposite(1,:),1/6*thirdorder_opposite(1,:),1/24*fourthorder_opposite(1,:)];
% dec_l_opposite=[LSTAR,firstorder_opposite(2,:),1/2*secondorder_opposite(2,:),1/6*thirdorder_opposite(2,:),1/24*fourthorder_opposite(2,:)];
% dec_c_opposite=[CSTAR,firstorder_opposite(3,:),1/2*secondorder_opposite(3,:),1/6*thirdorder_opposite(3,:),1/24*fourthorder_opposite(3,:)];


ssvals = [KSTAR,ZSTAR,ZSTAR,ZSTAR,ZSTAR,PSTAR];

rng(13466910+randomseq,'twister');
rho_zeta = normrnd(0,sigma_Z,[1 T]);
rho_zeta(1:5)=shock;

positiveshock = sigma_Z;
if shock == 88 %historical
    realTFP = readmatrix('TFPshocks.csv');
    rho_zeta(1:33+4) = realTFP(realTFP(:,1)>1980.5&realTFP(:,1)<2017.5,2);
    T = 37;
end
if shock == 55
    rho_zeta(1:5) = positiveshock;
    T = 33;
end

rng(123140+randomseq,'twister');
rho_P = normrnd(0,sigma_P,[1 T]);

rho_zeta(44:49) = 0;
rho_P(44:49) = 0;
        
k_sim=zeros([1 T]) + KSTAR*k0_mult;
Z_sim=zeros([1 T]) + ZSTAR^LAMBDAZ*exp(rho_zeta(1));
P_sim = zeros([1 T]) + P_func_greater_than_1(G,PSTAR,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
if startopposite
    P_sim = zeros([1 T]) + P_func_greater_than_1(G,Popposite,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
end

state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1)],T,1);
c_sim=zeros([1 T]) + decision_func_to_use(dec_c,state_vars(1,:),ssvals,sigma_Z,sigma_P);
l_sim=zeros([1 T]) + decision_func_to_use(dec_l,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% riskless_r_sim(1:T) = riskless_r_model2(state_vars(1,:),c_sim(1),l_sim(1),dec_c,dec_l,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym);
% riskless_r_sim_alt(1:T) = decision_func_to_use(dec_riskless_r,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% stock_sim(1:T) = decision_func_to_use(dec_stock,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% risky_r_sim(1:T) = decision_func_to_use(dec_risky_r,state_vars(1,:),ssvals,sigma_Z,sigma_P);
sophisticated_riskless_r = zeros([1 T]);
ktplus1 = zeros([1 T]);
eZtplus1 = zeros([1 T]);
ePtplus1 = zeros([1 T]);
e_ct_plus1 = zeros([1 T]);
e_cgrowth = zeros([1 T]);
euler = zeros([1 T]);

for i=2:T
    Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_zeta(i));

    if regimechanges&&((i>50&&i<=100)||(i>150&&i<=200))
        k_sim(i)=decision_func_to_use(dec_k_opposite,state_vars(i-1,:),ssvals,sigma_Z,sigma_P);

        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAPopposite,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
        
        state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),P_sim(i)];

        c_sim(i)=decision_func_to_use(dec_c_opposite,state_vars(i,:),ssvals,sigma_Z,sigma_P);
        l_sim(i)=decision_func_to_use(dec_l_opposite,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    
    else
        k_sim(i)=decision_func_to_use(dec_k,state_vars(i-1,:),ssvals,sigma_Z,sigma_P);

        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
        if i == 4 && shock == 88
            P_sim(i) = Popposite;
        end
        if i == 5 && shock == 88
            k_sim(i) = KSTAR*k0_mult;
        end

        state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),P_sim(i)];

        c_sim(i)=decision_func_to_use(dec_c,state_vars(i,:),ssvals,sigma_Z,sigma_P);
        l_sim(i)=decision_func_to_use(dec_l,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    end
    
%     riskless_r_sim_alt(i) = decision_func_to_use(dec_riskless_r,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     stock_sim(i) = decision_func_to_use(dec_stock,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     risky_r_sim(i) = decision_func_to_use(dec_risky_r,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     riskless_r_sim(i) = riskless_r_model2(state_vars(i,:),c_sim(i),l_sim(i),dec_c,dec_l,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym);
    ktplus1(i) = decision_func_to_use(dec_k,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    eZtplus1(i) = Z_sim(i)^LAMBDAZ;
    ePtplus1(i) = P_func_greater_than_1(G,P_sim(max(i,1)),LAMBDAP,eZtplus1(i),Z_sim(max(i,1)),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),MU,0);
    e_ct_plus1(i) = decision_func_to_use(dec_c,[ktplus1(i),Z_sim(i)^LAMBDAZ,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),ePtplus1(i)],ssvals,sigma_Z,sigma_P);
    [sophisticated_riskless_r(i),euler(i),e_cgrowth(i)] = riskless_r_model2([ktplus1(i),state_vars(i,2:6)],c_sim(i),l_sim(i),dec_c,dec_l,decision_func_to_use,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,ALFA,BETTA,G,GAMA,DELTA,ETA,SIGM,intertemporal_euler_sym);
end

premium = (G^SIGM/BETTA-e_cgrowth.*(1+sophisticated_riskless_r));
fprintf('Risk premium? %.10g\n',mean(premium))

w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA);
r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA);
y_sim = y_func(k_sim,l_sim,Z_sim,ALFA);
g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
profits_sim = ((P_sim-1)./P_sim).*y_sim;
riskless_r_naive = 1./(BETTA*(G*e_ct_plus1./c_sim).^(-SIGM))-1;

fprintf('median return %.10g\nmedian riskless rate %.10g\n',median(r_sim(2:T)),median(sophisticated_riskless_r(2:T)))
fprintf('mean return %.10g\nmean riskless rate %.10g\n',mean(r_sim(2:T)),mean(sophisticated_riskless_r(2:T)))

rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',y_sim',P_sim',Z_sim'];
if shock == 88
    rgwkcl_mat = rgwkcl_mat(5:37,:);
end
