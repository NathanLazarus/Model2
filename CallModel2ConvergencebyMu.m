DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
% LAMBDAP = 0.95;
LAMBDAZ = 0.909;
sigma_Z = 0.0109;
sigma_P = 0.02;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU = 0.37;
T=1000;
shock = 0.0127;


addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
import casadi.*


orders = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4];
LAMBDAPs = [0.8,0.8,0.8,0.8,0.95,0.95,0.95,0.95,0.8,0.8,0.8,0.8,0.95,0.95,0.95,0.95];
randomseqs = [2,2,2,2,2,2,2,2,13,13,13,13,13,13,13,13];
k0_mult=1;
startopposite = 0;
MultiplicativeU = 0;
    
N = length(orders);
output_vars=9;
output = NaN([T output_vars N]);
tic
parfor i = 1:N
    LAMBDAP = LAMBDAPs(i)
    order = orders(i)
    randomseq = randomseqs(i)
    [output(:,:,i)] = model2_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,randomseq,order);
end
toc
filename = "C:/Users/Nathan/Downloads/PerturbationMethods/Model_Outputs_May21.xlsx";
sheetloc = 'A6';
writematrix(output,filename,'Sheet',3,'Range',sheetloc)
writematrix(output(:,4,:),filename,'Sheet',1,'Range','A4')
writematrix(output(:,8,:),filename,'Sheet',2,'Range','A4')