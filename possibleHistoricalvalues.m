tic
N = 20000;
output = NaN([33 9 N]);
parfor j=1:N
    rng(123140+j,"twister");
    rho_P = normrnd(0,sigma_P,[1 T]);

    
    
    
    
    k_sim = zeros([1 T]) + KSTAR*k0_mult;
    Z_sim = zeros([1 T]) + ZSTAR^LAMBDAZ*exp(rho_Z(1));
    P_sim = zeros([1 T]) + P_func_greater_than_1(G,PSTAR,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
    if string(shock) == "historical_endogenous_P" && startopposite
        P_sim = zeros([1 T]) + P_func_greater_than_1(G,Popposite,LAMBDAP,Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),MU,rho_P(1));
    end

    state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1)],T,1);
    c_sim = zeros([1 T]) + decision_func_to_use(dec_c,state_vars(1,:),ssvals,sigma_Z,sigma_P);
    l_sim = zeros([1 T]) + decision_func_to_use(dec_l,state_vars(1,:),ssvals,sigma_Z,sigma_P);


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
            if i == T-(hardcode_irf_T-2) && string(shock_character_vector(1:3)) == "irf" && string(shock_character_vector(5)) == "K"
                k_sim(i) = k_sim(i-1)*0.5;
            end
            if i == 5 && (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
                k_sim(i) = KSTAR*k0_mult;
            end

            P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));

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

    end

    w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA);
    r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA);
    y_sim = y_func(k_sim,l_sim,Z_sim,ALFA);
    g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
    profits_sim = ((P_sim-1)./P_sim).*y_sim;

    profits = (P_sim-1)/P_sim*y_sim;

    rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',y_sim',P_sim',Z_sim'];
    if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
        rgwkcl_mat = rgwkcl_mat(5:37,:);
    end


    if string(shock_character_vector(1:3)) == "irf"
        rgwkcl_mat = rgwkcl_mat((T-(hardcode_irf_T-1)):T,:);
    end

    
    
    
    
    output(:,:,j) = rgwkcl_mat;
end
toc