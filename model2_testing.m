% model2_run.M
% Calls: model2.m num_eval.m  model2_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [rgwkcl_mat] = model2_testing(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,...
    sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,...
    startopposite,regimechanges,regime_change_frequency,randomseq,order,useopposite,shiftstart,startpoint)

hardcode_irf_T = 100;
eds_points = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/eds_points_for_irf.csv");
% clear all

use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients = false;

defaults = {0.08,0.32,0.98,1.014,...
    0.9,0.95,0.92,0.017,...
    0.03,0.086,0.5,0.3,...
    200,"none",1,1,...
    0,0,50,2,4,0,"none",1};

var = ["DELTA","ALFA","BETTA","G",...
    "SIGM","LAMBDAP","LAMBDAZ","sigma_Z",...
    "sigma_P","MU","FRISCHELAS","STEADYSTATEL",...
    "T","shock","k0_mult","MultiplicativeU",...
    "startopposite","regimechanges","regime_change_frequency","randomseq","order","useopposite","shiftstart","startpoint"];

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename));
folder
which(mfilename)
addpath(join([folder,'\','MyHelperFunctions'],""))
addpath(join([folder,'\','SchmittGroheUribeHelperFunctions'],""))


for i = 1:length(defaults)
    if ~exist(var{i},"var")
        if class(defaults{i}) == "string"
            eval(sprintf('%s = "%s";',var(i),defaults{i}))
        else
            eval(sprintf("%s = %g;",var(i),defaults{i}))
        end
    end
end

if MultiplicativeU
    u = multiplicative_u;
    utility_function_str = "Multiplicative";
else
    u = additive_u;
    utility_function_str = "Additive";
end

max_order_I_can_handle = 4;
nstates_plus_shocks = 8;
if use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients
    SchmittGrohe_decision_func = "SchmittGrohe";
    max_order_I_can_handle = 2;
else
    SchmittGrohe_decision_func = "";
end

decision_func_to_use = str2func(join(["decision_func_","order",max_order_I_can_handle,SchmittGrohe_decision_func],""));

LAMBDAPhigh = 0.95;
LAMBDAPlow = 0.8;
if LAMBDAP == LAMBDAPhigh
    LAMBDAPopposite = LAMBDAPlow;
end
if LAMBDAP == LAMBDAPlow
    LAMBDAPopposite = LAMBDAPhigh;
end

eta     = [0 0; 1 0; 0 0; 0 0; 0 0; 0 sigma_P/sigma_Z]; %Matrix defining driving force

ZSTAR = 1; %steady-state value of technology shock 
PSTAR = G^(1/(1-LAMBDAP)); %steady state markup

sym_labor_supply = laborsupply(u);
intertemporal_euler_ss = dupdcp_over_dudc(u,1);
intertemporal_euler_sym = dupdcp_over_dudc(u,0);

addpath("C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1")
import casadi.*

value_of_P_where_LSTAR_equals_STEADYSTATEL = G^(1/(1-LAMBDAPlow));

[~,~,~,~,~,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u);
% GAMA = 28.1677381532213;
[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model2_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTAR,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
% KSTAR = 0.835777225015446;
% LSTAR = 0.273094921033578;
% CSTAR = 0.312063920996547;
% RSTAR = little_r(KSTAR,LSTAR,PSTAR,1,ALFA,DELTA)
fprintf("{%.15g, %.15g, %.15g, %.15g, %.15g, %d}\n",KSTAR,CSTAR,LSTAR,GAMA,ETA,MultiplicativeU)
% [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,PSTAR,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u)
if startopposite
    if string(shock) == "historical_endogenous_P"
        Popposite = G^(1/(1-LAMBDAPopposite));
    end
    KSTARopposite = model2_ss_numeric(KSTAR,CSTAR,LSTAR,DELTA,ALFA,BETTA,G,Popposite,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
    k0_mult = KSTARopposite/KSTAR;
end



if use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients

    get_coefs_SchmittGrohe_Uribe
    
else

    for order_to_check_for_coefs=1:max_order_I_can_handle
        if order >= order_to_check_for_coefs
            eval(sprintf('order%1$d = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/MathematicaCoefs/",utility_function_str,"UtilityLambdaP",string(round(LAMBDAP*100)),"coefs",%1$d,".csv"],""));',order_to_check_for_coefs))
            eval(sprintf('order%1$d_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/MathematicaCoefs/",utility_function_str,"UtilityLambdaP",string(round(LAMBDAPopposite*100)),"coefs%1$d.csv"],""));',order_to_check_for_coefs))
        else
            eval(sprintf('order%1$d = zeros([3 nstates_plus_shocks^%1$d]);',order_to_check_for_coefs))
            eval(sprintf('order%1$d_opposite = zeros([3 nstates_plus_shocks^%1$d]);',order_to_check_for_coefs))
        end
    end

    dec_k=[KSTAR,order1(1,:),1/2*order2(1,:),1/6*order3(1,:),1/24*order4(1,:)];
    dec_l=[LSTAR,order1(2,:),1/2*order2(2,:),1/6*order3(2,:),1/24*order4(2,:)];
    dec_c=[CSTAR,order1(3,:),1/2*order2(3,:),1/6*order3(3,:),1/24*order4(3,:)];

    dec_k_opposite=[KSTAR,order1_opposite(1,:),1/2*order2_opposite(1,:),1/6*order3_opposite(1,:),1/24*order4_opposite(1,:)];
    dec_l_opposite=[LSTAR,order1_opposite(2,:),1/2*order2_opposite(2,:),1/6*order3_opposite(2,:),1/24*order4_opposite(2,:)];
    dec_c_opposite=[CSTAR,order1_opposite(3,:),1/2*order2_opposite(3,:),1/6*order3_opposite(3,:),1/24*order4_opposite(3,:)];
    
    if useopposite
        dec_k = dec_k_opposite;
        dec_l = dec_l_opposite;
        dec_c = dec_c_opposite;
    end
    
end

ssvals = [KSTAR,ZSTAR,ZSTAR,ZSTAR,ZSTAR,PSTAR];

rng(13466910+randomseq,"twister");
rho_Z = normrnd(0,sigma_Z,[1 T]);
if ~(string(shock) == "none"||string(shock) == "historical"|| string(shock) == "historical_endogenous_P")
    rho_Z(1:5)=shock;
end

rng(123140+randomseq,"twister");
rho_P = normrnd(0,sigma_P,[1 T]);


if (string(shock) == "historical" || string(shock) == "historical_endogenous_P") %historical
    T = 37;
    realTFP = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks.csv");
    rho_Z = zeros([1 T]) + realTFP(realTFP(:,1)>1980.5&realTFP(:,1)<2017.5,2);
    realProfits = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/ProfitShare.csv");
    Loop_P_Path = (1./(1-realProfits(realProfits(:,1)>1980.5&realProfits(:,1)<2017.5,2)))';
    simulated_P_implied_by_historical_Z = Loop_P_Path;
    historical_Z_path = zeros([1 T]) + ZSTAR^LAMBDAZ*exp(rho_Z(1));
    rho_P = zeros([1 T]) + log(Loop_P_Path(1)/P_func(G,PSTAR,LAMBDAP,historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),MU,0));
    for i=2:T
        historical_Z_path(i)=historical_Z_path(i-1)^LAMBDAZ*exp(rho_Z(i));
        historical_P_implied_by_Z = P_func(G,Loop_P_Path(max(i-1,1)),LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
        simulated_P_implied_by_historical_Z(i) = P_func(G,simulated_P_implied_by_historical_Z(max(i-1,1)),LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
        rho_P(i) = log(Loop_P_Path(i)/historical_P_implied_by_Z);
    end
end

shock_character_vector = char(shock);

if string(shock_character_vector(1:3)) == "irf"
    rho_Z = zeros([1 T]);
    rho_P = zeros([1 T]);
    if string(shock_character_vector(5)) == "Z"
        rho_Z(T-(hardcode_irf_T-2)) = sigma_Z;
    end
    if string(shock_character_vector(5)) == "z"
        rho_Z(T-(hardcode_irf_T-2)) = sigma_Z;
        rho_Z(T-(hardcode_irf_T-3)) = log(1/exp(sigma_Z)^LAMBDAZ);
    end
    if string(shock_character_vector(5)) == "P"
        rho_P(T-(hardcode_irf_T-2)) = sigma_P;
    end
    if string(shock_character_vector(5)) == "K"
        k1_mult = 0.5;
    end
end
        

%Initializing Arrays

k_sim = zeros([1 T]) + KSTAR*k0_mult;
Z_sim = zeros([1 T]) + ZSTAR^LAMBDAZ*exp(rho_Z(1));
P_sim = zeros([1 T]) + P_func_greater_than_1(G,PSTAR,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
if string(shock) == "historical_endogenous_P" && startopposite
    P_sim = zeros([1 T]) + P_func_greater_than_1(G,Popposite,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
end

state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1)],T,1);
c_sim = zeros([1 T]) + decision_func_to_use(dec_c,state_vars(1,:),ssvals,sigma_Z,sigma_P);
l_sim = zeros([1 T]) + decision_func_to_use(dec_l,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% riskless_r_sim(1:T) = riskless_r_model2(state_vars(1,:),c_sim(1),l_sim(1),dec_c,dec_l,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym);
% riskless_r_sim_alt(1:T) = decision_func_to_use(dec_riskless_r,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% stock_sim = ([1 T]) + decision_func_to_use(dec_stock,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% risky_r_sim = ([1 T]) + decision_func_to_use(dec_risky_r,state_vars(1,:),ssvals,sigma_Z,sigma_P);
% sophisticated_riskless_r = zeros([1 T]);
% ktplus1 = zeros([1 T]);
% eZtplus1 = zeros([1 T]);
% ePtplus1 = zeros([1 T]);
% e_ct_plus1 = zeros([1 T]);
% e_cgrowth = zeros([1 T]);
% euler = zeros([1 T]);

for i=2:T
    Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_Z(i));

    if regimechanges && mod(floor((i-1)/regime_change_frequency),2) == 1
        k_sim(i)=decision_func_to_use(dec_k_opposite,state_vars(i-1,:),ssvals,sigma_Z,sigma_P);

        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAPopposite,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
        
        state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),P_sim(i)];

        c_sim(i)=decision_func_to_use(dec_c_opposite,state_vars(i,:),ssvals,sigma_Z,sigma_P);
        l_sim(i)=decision_func_to_use(dec_l_opposite,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    
    else
        k_sim(i)=decision_func_to_use(dec_k,state_vars(i-1,:),ssvals,sigma_Z,sigma_P);
        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));

        if i == T-(hardcode_irf_T-2) && string(shock_character_vector(1:3)) == "irf" && string(shock_character_vector(5)) == "K"
            k_sim(i) = k_sim(i-1)*0.5;
        end
        if i == T-(hardcode_irf_T-1) && string(shock_character_vector(1:3)) == "irf" && contains(shiftstart,"k")
            k_sim(i) = eds_points(startpoint,1);
        end
        if i < T-(hardcode_irf_T-1) && string(shock_character_vector(1:3)) == "irf" && contains(shiftstart,"P")
            P_sim(i) = eds_points(startpoint,2);
        end
        if i < T-(hardcode_irf_T-1) && string(shock_character_vector(1:3)) == "irf" && contains(shiftstart,"Z")
            Z_sim(max(i-3,1)) = eds_points(startpoint,6);
            Z_sim(max(i-2,1)) = eds_points(startpoint,5);
            Z_sim(max(i-1,1)) = eds_points(startpoint,4);
            Z_sim(i) = eds_points(startpoint,3);
        end
        if i == 5 && (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
            k_sim(i) = KSTAR*k0_mult;
        end
        
        if i == 4 && string(shock) == "historical_endogenous_P" && startopposite
            P_sim(i) = Popposite;
        end
        if i == 4 && string(shock) == "historical_endogenous_P" && ~startopposite
            P_sim(i) = PSTAR;
        end

        state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),P_sim(i)];

        c_sim(i)=decision_func_to_use(dec_c,state_vars(i,:),ssvals,sigma_Z,sigma_P);
        l_sim(i)=decision_func_to_use(dec_l,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    end
    
%     riskless_r_sim_alt(i) = decision_func_to_use(dec_riskless_r,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     stock_sim(i) = decision_func_to_use(dec_stock,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     risky_r_sim(i) = decision_func_to_use(dec_risky_r,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     riskless_r_sim(i) = riskless_r_model2(state_vars(i,:),c_sim(i),l_sim(i),dec_c,dec_l,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym);
%     ktplus1(i) = decision_func_to_use(dec_k,state_vars(i,:),ssvals,sigma_Z,sigma_P);
%     eZtplus1(i) = Z_sim(i)^LAMBDAZ;
%     ePtplus1(i) = P_func_greater_than_1(G,P_sim(max(i,1)),LAMBDAP,eZtplus1(i),Z_sim(max(i,1)),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),MU,0);
%     e_ct_plus1(i) = decision_func_to_use(dec_c,[ktplus1(i),Z_sim(i)^LAMBDAZ,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),ePtplus1(i)],ssvals,sigma_Z,sigma_P);
%     [sophisticated_riskless_r(i),euler(i),e_cgrowth(i)] = riskless_r_model2([ktplus1(i),state_vars(i,2:6)],c_sim(i),l_sim(i),dec_c,dec_l,decision_func_to_use,LAMBDAZ,LAMBDAP,MU,ssvals,sigma_Z,sigma_P,ALFA,BETTA,G,GAMA,DELTA,ETA,SIGM,intertemporal_euler_sym);
end

w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA);
r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA);
y_sim = y_func(k_sim,l_sim,Z_sim,ALFA);
g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
profits_sim = ((P_sim-1)./P_sim).*y_sim;
% riskless_r_naive = 1./(BETTA*(G*e_ct_plus1./c_sim).^(-SIGM))-1;

profits = (P_sim-1)/P_sim*y_sim;
% stockreturns = (stock_sim(2:T)+profits(2:T))./stock_sim(1:(T-1));
% premium = (G^SIGM/BETTA-e_cgrowth.*(1+sophisticated_riskless_r));
% fprintf("Risk premium? %.10g\n",mean(premium))
% 
% fprintf("%.10g\n%.10g\n%.10g\n",median(r_sim(2:T)),median(riskless_r_naive(2:T)),median(sophisticated_riskless_r(2:T)))
% fprintf("%.10g\n%.10g\n%.10g\n",mean(r_sim(2:T)),mean(riskless_r_naive(2:T)),mean(sophisticated_riskless_r(2:T)))

rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',y_sim',P_sim',Z_sim'];
if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
    rgwkcl_mat = rgwkcl_mat(5:37,:);
end


if string(shock_character_vector(1:3)) == "irf"
    rgwkcl_mat = rgwkcl_mat((T-(hardcode_irf_T-1)):T,:);
end

% 
% filename = join(["C:/Users/Nathan/Downloads/PerturbationMethods/IRFs.xlsx"],"");
% sheetloc = 'A5'; %B4
% writematrix(rgwkcl_mat,filename,'Sheet',1,'Range',sheetloc)
