% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [rgwkcl_mat] = model2_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU)

defaults = [0.12,0.32,0.96,1.018,0.9,0.95,0.896,0.0072,0.005,0.35,0.5,0.3,10,0.015,1,0];
var={'DELTA','ALFA','BETTA','G','SIGM','LAMBDAP','LAMBDAZ','sigma_Z','sigma_P','MU','FRISCHELAS','STEADYSTATEL','T','shock','k0_mult','MultiplicativeU'};

for i = 1:length(defaults)
    if ~exist(var{i},'var')
        eval(sprintf('%s = %g',var{i},defaults(i)))
    end
end

if MultiplicativeU
    u = multiplicative_u;
else
    u = additive_u;
end

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model2(u);

eta     = [0 0; 1 0; 0 0; 0 0; 0 0; 0 0; 0 sigma_P/sigma_Z]; %Matrix defining driving force
T_Psims = T;

% impulse response functions setup
irf=0;
Z_shock=sigma_Z;
T_irf = T;

% simulations setup
simulations=1;
irf = 0;

ZSTAR = 1; %steady-state value of technology shock 
PSTAR = G^(1/(1-LAMBDAP)); %steady state markup

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,PSTAR,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,u);


k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR; Z1 = ZSTAR; Z2 = ZSTAR; Z3 = ZSTAR; Z4 = ZSTAR; P = PSTAR;
kp=k; cp=c; lp=l; Zp=Z; Z1p = Z; Z2p = Z; Z3p = Z; Z4p = Z; Pp = P;

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
    dec_c=[LSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];

else
    dec_k=[KSTAR,hx(1,:),zeros([1 nstate^2+1])]; 
    dec_l=[LSTAR,gx(1,:),zeros([1 nstate^2+1])];
    dec_c=[CSTAR,gx(2,:),zeros([1 nstate^2+1])];
    
end

ssvals = [KSTAR,ZSTAR,ZSTAR,ZSTAR,ZSTAR,ZSTAR,PSTAR];
if irf
    k(1)=KSTAR;
    c(1)=CSTAR;
    Z(1)=ZSTAR;
    l(1)=LSTAR;
    P(1)=PSTAR;
    
    state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1)],T_irf,1);

    for i=2:40
        k(i)=decision_func(dec_k,state_vars(i-1,:),ssvals,sigma_Z);
        if i==2
            Z(i) = Z(i-1)^LAMBDAZ*exp(Z_shock);
        else
            Z(i) = Z(i-1)^LAMBDAZ;
        end
        P(i) = P_func(G,P(max(i-1,1)),LAMBDAP,Z(max(i-1,1)),Z(max(i-2,1)),Z(max(i-3,1)),Z(max(i-4,1)),Z(max(i-5,1)),MU,rho_P(i));        
        
        state_vars(i,:) = [k(i),Z(i),Z(max(i-1,1)),Z(max(i-2,1)),Z(max(i-3,1)),Z(max(i-4,1)),P(i)];
        c(i)=decision_func(dec_c,state_vars(i,:),ssvals,sigma_Z);
        l(i)=decision_func(dec_c,state_vars(i,:),ssvals,sigma_Z);
    end

    w = w_func(k,l,P,Z,ALFA);
    r = little_r(k,l,P,Z,ALFA,DELTA);

    % figure(1)
    % subplot(2,3,1)
    % plot(k)
    % title('Capital ($k$)','Interpreter','latex')
    % subplot(2,3,2)
    % plot(c)
    % title('Consumption ($c$)','Interpreter','latex')
    % subplot(2,3,3)
    % plot(l)
    % title('Labor ($l$)','Interpreter','latex')
    % subplot(2,3,4)
    % plot(Z)
    % title('Productivity shock ($\zeta$)','Interpreter','latex')
    % subplot(2,3,5)
    % plot(r)
    % title('Interest rate ($r$)','Interpreter','latex')
    % subplot(2,3,6)
    % plot(w)
    % title('Wage rate ($w$)','Interpreter','latex')
    % print(['IRF_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    %close(1)
end

rng(13466910,'twister');
rho_zeta = normrnd(0,sigma_Z,[1 max(T,T_Psims)]);
rho_zeta(1:5)=shock;

rng(20,'twister');
rho_P = normrnd(0,sigma_P,[1 T_Psims]);
        
if simulations
    %fprintf('\n mean(rho_zeta)=%g, estimated=%g\n',0, mean(rho_zeta))
    %fprintf('sigma_Z=%g, estimated=%g\n\n', sigma_Z, std(rho_zeta))
    
    % Start from the non-stochastic steady state
    k_sim(1:T)=KSTAR*k0_mult;
    Z_sim(1:T)=ZSTAR*exp(rho_zeta(1));
    P_sim(1:T) = P_func(G,PSTAR,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
    
    state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1)],T,1);
    c_sim(1:T)=decision_func(dec_c,state_vars,ssvals,sigma_Z);
    l_sim(1:T)=decision_func(dec_l,state_vars,ssvals,sigma_Z);
    
    for i=2:T
        k_sim(i)=decision_func(dec_k,state_vars(i-1,:),ssvals,sigma_Z);
        Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_zeta(i));
        P_sim(i) = P_func(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),Z_sim(max(i-5,1)),MU,rho_P(i));        
        
        state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),P_sim(i)];
        c_sim(i)=decision_func(dec_c,state_vars(i,:),ssvals,sigma_Z);
        l_sim(i)=decision_func(dec_c,state_vars(i,:),ssvals,sigma_Z);
    end
    
    w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA);
    r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA);
    y_sim = y_func(k_sim,l_sim,Z_sim,ALFA);
    g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
    
    % figure(2)
    % subplot(2,3,1)
    % plot(k_sim)
    % title('Capital ($k$)','Interpreter','latex')
    % subplot(2,3,2)
    % plot(c_sim)
    % title('Consumption ($c$)','Interpreter','latex')
    % subplot(2,3,3)
    % plot(l_sim)
    % title('Labor ($l$)','Interpreter','latex')
    % subplot(2,3,4)
    % plot(Z_sim)
    % title('Productivity shock ($\zeta$)','Interpreter','latex')    
    % subplot(2,3,5)
    % plot(r_sim)
    % title('Interest rate ($r$)','Interpreter','latex')
    % subplot(2,3,6)
    % plot(w_sim)
    % title('Wage rate ($w$)','Interpreter','latex')
    % print(['SIM_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    %close(2)
end

rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',P_sim',Z_sim'];
if exist('filename','var')   
    writematrix(rgwkcl_mat,filename,'Sheet',1,'Range',sheetloc)
end