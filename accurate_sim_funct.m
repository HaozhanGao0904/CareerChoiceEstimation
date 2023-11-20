% this function is to simulate data from parameter
% simulation is more accurate since it uses accurate EMAX (no
% interpolation)
% simulation function

function sim_data = accurate_sim_funct(R_here,fp,param)

%%%%%%%%%%%% unpack parameters (40 in total) %%%%%%%%%%%%
r1 = param(1);
r2 = param(2);
e11 = param(3);
e12 = param(4);
e13 = param(5);
e14 = param(6);
e15 = param(7);
    % i dropped the quadratic term for bc exp for wc job
e21 = param(8);
e22 = param(9);
e23 = param(10);
    % i dropped the quadratic term for wc exp for bc job
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

%%%%%%%%%%%%%%%%% obtain EMAX matrix %%%%%%%%%%%%%%%%%%%%

EMAX = get_accurate_EMAX(fp, param);

%%%%%%%%%%%%%%%%% solve forward %%%%%%%%%%%%%%%%%%%%

% draw random shocks
mean_vec = [0 0 0 0];

cov_matrix = [s11,  s12,   0,   0;
              s12,    s22, 0,   0;
              0,    0,   s33*10^6, 0;
              0,    0,   0,   s44*10^6];

rng(614);
fp.forward_draws = zeros(fp.R_forward_max, fp.C, fp.T);
for page = 1:fp.T
    fp.forward_draws(:,:,page) = mvnrnd(mean_vec, cov_matrix, fp.R_forward_max);
end
   

sim_data = nan(fp.T*R_here, fp.data_cols);

index = 0;

plow = [1-plow2-plow3-plow4, plow2, plow3, plow4];
phigh = [1-phigh2-phigh3-phigh4, phigh2, phigh3, phigh4];

rng('default');
for r = 1:R_here
    % initial work exp and edu yrs
    x1 = 0;
    x2 = 0;
    g = fp.ini_edu_male(r); 
    j = fp.ini_edu_male(r)-6; % page number
    
    % realize her type

    if g <= 9
        type = find(rand <= cumsum(plow), 1);
    else
        type = find(rand <= cumsum(phigh), 1);
    end

    % realize her endowment
    if type == 1
        r1_h = r1;
        r2_h = r2;
        e30_h = e30;
        e40_h = e40;
    elseif type == 2
        r1_h = r1 - t12;
        r2_h = r2 - t22;
        e30_h = e30 - t32;
        e40_h = e40 - t42;
    elseif type == 3
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
    
    % indicator for job switching
    indc1 = 0;
    indc2 = 0;

    for t = 1:fp.T
        age = 16 + t;
        x1_index = x1 + 1;
        x2_index = x2 + 1;
        g_index = g + 1 - fp.ini_edu_male(r);
      
        % 1. white-collar
        x1_index_next_period = x1_index + 1;
        x2_index_next_period = x2_index;
        g_index_next_period = g_index;

        log_wage_wc = 0.01*(e11*g + e12*g^2 + e13*x1 + e14*x1^2 + e15*x2) + fp.forward_draws(r,1,t);

        u_white_collar = 1000*r1_h*exp(log_wage_wc) + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, j);

        % 2. blue-collar
        x1_index_next_period = x1_index;
        x2_index_next_period = x2_index + 1;
        g_index_next_period = g_index;

        log_wage_bc = 0.01*(e21*g + e22*g^2 + e23*x1 + e25*x2...
                   +e26*x2^2) + fp.forward_draws(r,2,t);
        u_blue_collar = 1000*r2_h*exp(log_wage_bc) + delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, j);



        % 3. school
        x1_index_next_period = x1_index;
        x2_index_next_period = x2_index;
        g_index_next_period = g_index + 1;

        g_indc1 = 0; % indicator of whether g>12
        g_indc2 = 0; % indicator of whether g>16
        g_indc1(g > 12) = 1;
        g_indc2(g > 16) = 1;

        u_school = 10000*(e30_h - c1*g_indc1 - c2*g_indc2) + fp.forward_draws(r,3,t) + ...
                               delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, j);

        % 4. home
        x1_index_next_period = x1_index;
        x2_index_next_period = x2_index;
        g_index_next_period = g_index;

        u_home = 10000*e40_h + fp.forward_draws(r,4,t) + ...
                 delta.*EMAX(x1_index_next_period, x2_index_next_period, g_index_next_period, t+1, j);

        % optimal choices:

        if g >= 20 || t >= 20
            if u_white_collar == max([u_white_collar, u_blue_collar, u_home])
                d_opt = 1;
                log_wage_obs = log_wage_wc;
                wage = 1000*r1_h*exp(log_wage_obs);
            elseif u_blue_collar == max([u_white_collar, u_blue_collar, u_home])
                d_opt = 2;
                log_wage_obs = log_wage_bc;
                wage = 1000*r2_h*exp(log_wage_obs);
            else
                d_opt = 4;
                wage = nan;
            end 
        else
            if u_white_collar == max([u_white_collar, u_blue_collar, u_school, u_home])
                d_opt = 1;
                log_wage_obs = log_wage_wc;
                wage = 1000*r1_h*exp(log_wage_obs);
            elseif u_blue_collar == max([u_white_collar, u_blue_collar, u_school, u_home])
                d_opt = 2;
                log_wage_obs = log_wage_bc;
                wage = 1000*r2_h*exp(log_wage_obs);
            elseif u_school == max([u_white_collar, u_blue_collar, u_school, u_home])
                d_opt = 3;
                wage = nan;
            else
                d_opt = 4;
                wage = nan;
            end
        end

        index = index + 1;


        sim_data(index,:) = [r t d_opt x1 x2 g log(wage)];

        % information updating
        if d_opt == 1
            x1 = x1 + 1;
            indc1 = 0;
            indc2 = 1;
        elseif d_opt == 2
            x2 = x2 + 1;
            indc1 = 1;
            indc2 = 0;
        elseif d_opt == 3
            g = g + 1;
            indc1 = 0;
            indc2 = 0;
        else
            indc1 = 1;
            indc2 = 1;
        end

    end % t loop
    

end % r loop
