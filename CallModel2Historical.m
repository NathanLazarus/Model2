clear all

DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
LAMBDAP = 0.8;
LAMBDAZ = 0.92;
sigma_Z = 0.017;
sigma_P = 0.03;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU = 0.086;
T=1000;
k0_mult=1;
startopposite = 0;
regimechanges = 0;
regime_change_frequency = 50;
order = 4;
randomseq = 2;


MultiplicativeU = 1;

addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
import casadi.*

if MultiplicativeU
    utility_function_str = "Multiplicative";
else
    utility_function_str = "Additive";
end

shocks = ["historical_postwar_trend","historical","historical_endogenous_P_postwar_trend","historical_endogenous_P"];

    
N = length(shocks);
output_vars=10;
output = NaN([33 output_vars N]);
tic
parfor i = 1:N
    shock = shocks(i);
    [output(:,:,i)] = model2_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,regimechanges,regime_change_frequency,randomseq,order);
end
toc

writematrix(output(:,:,1:2),"C:/Users/Nathan/Downloads/PerturbationMethods/Historical_Model_Outputs.xlsx",'Sheet',1,'Range','B4')
writematrix(output(:,:,3:4),"C:/Users/Nathan/Downloads/PerturbationMethods/Historical_Model_Outputs.xlsx",'Sheet',2,'Range','B4')
