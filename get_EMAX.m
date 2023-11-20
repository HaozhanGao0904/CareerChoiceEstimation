% this function will obtain the EMAX matrix

% 1. for selected spaces, I use simulation to obtain their "true" EMAX
% 2. using simulated EMAX,MAXE and V_m, run regression to obtain coefficients
% 3. for the rest of the spaces, I use interpolation to get their EMAX

function EMAX = get_EMAX(fp, param)

% unpack parameters (40 in total)
r1 = param(1);
r2 = param(2);
e11 = param(3);
e12 = param(4);
e13 = param(5);
e14 = param(6);
e15 = param(7);

e21 = param(8);
e22 = param(9);
e23 = param(10);

e25 = param(11);
e26 = param(12);
e30 = param(13);
c1 = param(14);
c2 = param(15);
e40 = param(16);
delta = param(17);

s11 = param(18);
s22 = param(19);
s33 = param(20);
s44 = param(21);
s12 = param(22);
%{
plow2 = param(23);
plow3 = param(24);
plow4 = param(25);
phigh2 = param(26);
phigh3 = param(27);
phigh4 = param(28);
%}
t12 = param(29);
t22 = param(30);
t32 = param(31);
t42 = param(32);
t13 = param(33);
t23 = param(34);
t33 = param(35);
t43 = param(36);
t14 = param(37);
t24 = param(38);
t34 = param(39);
t44 = param(40);


% define a EMAX matrix with placeholders of zeros
EMAX = zeros(fp.T+1,fp.T+1,fp.T+1,fp.T+1,fp.E,fp.Type);

% random shock draw
mean_vec = [0 0 0 0];

cov_matrix = [s11    s12      0        0;
              s12    s22    0        0;
              0      0      s33*1e6  0;
              0      0      0        s44*1e6];

rng('default');
fp.EMAX_draws = mvnrnd(mean_vec, cov_matrix, fp.R_EMAX);




for ini_edu_index = 1:fp.E % loop over initial education

    [space, selectedspace] = select_space(fp, ini_edu_index);

    for type_index = 1:fp.Type % loop over type space

    % define intercepts for heterogeneous types
        if type_index == 1
            r1_h = r1;
            r2_h = r2;
            e30_h = e30;
            e40_h = e40;
        elseif type_index == 2
            r1_h = r1 - t12;
            r2_h = r2 - t22;
            e30_h = e30 - t32;
            e40_h = e40 - t42;
        elseif type_index == 3
            r1_h = r1 - t13;
            r2_h = r2 - t23;
            e30_h = e30 - t33;
            e40_h = e40 - t43; 
        else
            r1_h = r1 - t14;
            r2_h = r2 - t24;
            e30_h = e30 - t34;
            e40_h = e40 - t44;
        end
    
    % compute EMAX for each period

    % use simulation to get EMAX for selected space
    % we also create lists to store EMAX, MAXE, Vm (m = 1,..,4) that are used for regression
    for t = fp.T:-1:2
        len = length(selectedspace{t});

        EMAX_list = nan(len,1);
        MAXE_list = nan(len,1);
        V1_list = nan(len,1);
        V2_list = nan(len,1);
        V3_list = nan(len,1);
        V4_list = nan(len,1);

        for i = 1:len % for every state at time t
            state = selectedspace{t}{i}; % get states
    
            % extract state variables
            x1 = state(1);
            x2 = state(2);
            g = state(3);
            
            % create indices
            x1_index = x1 + 1;
            x2_index = x2 + 1;
            g_index = g - 6 - ini_edu_index + 1; % since g = 6 + ini_edu_index + (g_index - 1)

            % creat placeholder lists
            
            u_max = nan(fp.R_EMAX, 1);
            u_wc_list = nan(fp.R_EMAX, 1);
            u_bc_list = nan(fp.R_EMAX, 1);
            u_school_list = nan(fp.R_EMAX, 1);
            u_home_list = nan(fp.R_EMAX, 1);

            for r = 1:fp.R_EMAX
                
         
                % 1. white-collar
                x1_index_next_period = x1_index + 1;
                x2_index_next_period = x2_index;
                g_index_next_period = g_index;
            
                log_wage_wc = 0.01*(e11*g + e12*g^2 + e13*x1 + e14*x1^2 + e15*x2) + fp.EMAX_draws(r, 1);
                u_white_collar = 1000*r1_h*exp(log_wage_wc) + ...
                                    + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
            
                 % 2. blue-collar
                 x1_index_next_period = x1_index;
                 x2_index_next_period = x2_index + 1;
                 g_index_next_period = g_index;
            
                 log_wage_bc = 0.01*(e21*g + e22*g^2 + e23*x1 + e25*x2...
                                    +e26*x2^2) + fp.EMAX_draws(r, 2);
                 u_blue_collar = 1000*r2_h*exp(log_wage_bc) + ...
                                    delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
            
                 % 3. school
                 x1_index_next_period = x1_index;
                 x2_index_next_period = x2_index;
                 g_index_next_period = g_index + 1;
            
                 g_indc1 = 0; % indicator of whether g>12
                 g_indc2 = 0; % indicator of whether g>16
                 g_indc1(g > 12) = 1;
                 g_indc2(g > 16) = 1;
            
                 u_school = 10000*(e30_h - c1*g_indc1 - c2*g_indc2) + fp.EMAX_draws(r, 3) + ...
                                delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
            
                 % 4. home
                 x1_index_next_period = x1_index;
                 x2_index_next_period = x2_index;
                 g_index_next_period = g_index;
            
                 u_home = 10000*e40_h + fp.EMAX_draws(r, 4) + ...
                            delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
            
                 % max utility among the 4 choices
                 if g > 20
                    u_max(r) = max([u_white_collar, u_blue_collar, u_home]);               
                 else
                    u_max(r) = max([u_white_collar, u_blue_collar, u_school, u_home]);
                 end

                 % store utility for 4 choices
                 u_wc_list(r) = u_white_collar;
                 u_bc_list(r) = u_blue_collar;
                 u_school_list(r) = u_school;
                 u_home_list(r) = u_home;
                 
            end % end of r loop

            % fill in EMAX matrix for selected space
            EMAX(x1_index, x2_index, g_index, t, ini_edu_index, type_index) = mean(u_max); 
            
            % store values in lists that will be used for regression
            EMAX_list(i) = mean(u_max);
            MAXE_list(i) = max([ mean(u_wc_list), mean(u_bc_list), mean(u_school_list), mean(u_home_list) ]);
            V1_list(i) = mean(u_wc_list);
            V2_list(i) = mean(u_bc_list);
            V3_list(i) = mean(u_school_list);
            V4_list(i) = mean(u_home_list);

        end % end of i loop
        
        % now we need to fill in the rest states in EMAX matrix
        spacelen = length(space{t});

        if len < spacelen % if ||selectedspace|| < ||space||  
                % run regression to get coefficients for interpolation
                % coeffients are time-specific
                % EMAX-MAXE = pi_0 + \sum_{m=1}^4 pi_1m (MAXE-V_m) + ...
                % \sum_{m=1}^4 pi_2m(MAXE-Vm)^{1/2}
                y = EMAX_list - MAXE_list;
                X = [ones(len,1), MAXE_list-V1_list, MAXE_list-V2_list, MAXE_list-V3_list, MAXE_list-V4_list, ...
                    (MAXE_list-V1_list).^0.5, (MAXE_list-V2_list).^0.5, (MAXE_list-V3_list).^0.5, (MAXE_list-V4_list).^0.5];
                pi_list = regress(y,X);
                %disp(pi_list)
               
               % now use coefficients to interpolate total state space of t
        
               for j = 1:spacelen
                    state = space{t}{i}; % get states from the whole space

                    % extract state variables
                    x1 = state(1);
                    x2 = state(2);
                    g = state(3); 
                    
                    % create indices
                    x1_index = x1 + 1;
                    x2_index = x2 + 1;
                    g_index = g - 6 - ini_edu_index + 1; % since g = 6 + ini_edu_index + (g_index - 1)
                      
                    if EMAX(x1_index, x2_index, g_index, t, ini_edu_index, type_index) == 0 % for all the unfilled states
                    % compute indepent variables: MAXE, V_m directly
                    % since we assume shocks have mean 0    
    
                    % 1: white-collar
                    x1_index_next_period = x1_index + 1;
                    x2_index_next_period = x2_index;
                    g_index_next_period = g_index;
    
                    log_wage_wc = 0.01*(e11*g + e12*g^2 + e13*x1 + e14*x1^2 + e15*x2);
                    V1 = 1000*r1_h*exp(log_wage_wc) + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
                    
                    % 2. blue-collar
                    x1_index_next_period = x1_index;
                    x2_index_next_period = x2_index + 1;
                    g_index_next_period = g_index;
                
                    log_wage_bc = 0.01*(e21*g + e22*g^2 + e23*x1 + e25*x2 + e26*x2^2);
                    V2 = 1000*r2_h*exp(log_wage_bc) + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
    
                    % 3. school    
    
                    x1_index_next_period = x1_index;
                    x2_index_next_period = x2_index;
                    g_index_next_period = g_index + 1;
                
                    g_indc1 = 0; % indicator of whether g>12
                    g_indc2 = 0; % indicator of whether g>16
                    g_indc1(g > 12) = 1;
                    g_indc2(g > 16) = 1;
                
                    V3 = 10000*(e30_h - c1*g_indc1 - c2*g_indc2) + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
                
                     % 4. home
                    x1_index_next_period = x1_index;
                    x2_index_next_period = x2_index;
                    g_index_next_period = g_index;
                
                    V4 = 10000*e40_h + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, ini_edu_index, type_index);
                    
    
                    MAXE = max([V1,V2,V3,V4]);

                    % compute EMAX using coefficients
                    EMAX(x1_index, x2_index, g_index, t, ini_edu_index, type_index) = ...
                        MAXE + pi_list(1) + pi_list(2)*(MAXE-V1) + pi_list(3)*(MAXE-V2) + pi_list(4)*(MAXE-V3) + pi_list(5)*(MAXE-V4)+...
                        pi_list(6)*((MAXE-V1)^0.5) + pi_list(7)*((MAXE-V2)^0.5) + pi_list(8)*((MAXE-V3)^0.5) +pi_list(9)*((MAXE-V4)^0.5);
        
                    end % end of if condition
               end
        end
    end % end of t loop
    end % end of type loop
end % end of initial education loop