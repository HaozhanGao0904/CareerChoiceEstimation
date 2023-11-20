% this function will generate EMAX matrix without using interpolation
% so it takes a longer time but the data is less noisy
function EMAX = get_accurate_EMAX(fp, param)

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
plow2 = param(23);
plow3 = param(24);
plow4 = param(25);
phigh2 = param(26);
phigh3 = param(27);
phigh4 = param(28);
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

cov_matrix = [s11    0      0        0;
              0      s22    0        0;
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
    % use simulation to get EMAX for entire space
    for t = fp.T:-1:2
        len = length(space{t});

        for i = 1:len % for every state at time t
            state = space{t}{i}; % get states
    
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

        end % end of i loop
    end % end of t loop
    end % end of type loop
end % end of initial education loop