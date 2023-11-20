%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Career Choice of the New Generation             %
%               Haozhan Gao(UW-Madison), 2023               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off', 'all'); % suppress the warning for regressions in the moment function

%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%

fp.T = 50; % total time period, from age 16 to 65
fp.D = 11; % recorded time period in real data (used for graph function)
fp.C = 4; % number of available choices each period
fp.N = 1000; % number of agents of simulation

% probability of different initial education level (obtained from data)
% male:
fp.M7 = 0.0079; % propotion to have initial education of 7 years for male
fp.M8 = 0.0383;
fp.M9 = 0.3261;
fp.M10 = 0.5059;
fp.M11 = 0.1218;

% female:
fp.F7 = 0.0026; % propotion to have initial education of 7 years for female
fp.F8 = 0.0313;
fp.F9 = 0.2889;
fp.F10 = 0.5265;
fp.F11 = 0.1506;


% estimation/simulation parameters

fp.J = 40; % number of parameters
fp.E = 5; % 5 possible initial education levels
fp.Type = 4; % 4 types
fp.M = 12024; % number of moments
fp.R_EMAX = 100; % number of simulation agents to use their mean to compute EMAX
fp.Interp = 1000; % number of points that are selected for interpolation
fp.R_forward_max = fp.N; % max forward draws
fp.B = 30; % bootstrap reps for calculating covariance
%fp.W = eye(fp.M); % weighting matrix


% data column indices
fp.col_r = 1; % 1st column donotes observation
fp.col_t = 2; % 2nd column denotes time period
fp.col_d = 3; % 3rd column denotes decision (from 1 to 4)
              % 1: white-collar; 2: blue-collar; 3: school; 4: home
fp.col_x1 = 4; % 4th column denotes white-collar work exp
fp.col_x2 = 5; % 5th column denotes blue-collar work exp
fp.col_g = 6; % 6th column denotes education years
fp.col_lnw = 7; % 7th column denotes lnwage

fp.data_cols = 7; % total columns in the table

%%%%%%%%%%%%%%%%%% agent parameters %%%%%%%%%%%%%%%%%%
%
% conditional type probabilities
% for ini_edu <= 9 (lower eduction)
plow2 = 0.2396; % probability of being type 2 when education is low
plow3 = 0.5015;
plow4 = 0.0838;

% for ini_edu >= 10 (higher education)
phigh2 = 0.4409; % probability of being type 2 when education is high
phigh3 = 0.4876;
phigh4 = 0.0329;

% parameters for labor rental price
r1_true = 9.4; % white-collar intercept for type 1
t12 = 0.0668; % deviation of type 2 from type 1 in white-collar job
t13 = 0.2221; % deviation of type 3 from type 1 in white-collar job
t14 = 0.2998; % deviation of type 4 from type 1 in white-collar job

r2_true = 10.6; % blue-collar intercept for type 1
t22 = 0.2996; % deviation of type 2 from type 1 in blue-collar job
t23 = -0.1223; % deviation of type 3 from type 1 in blue-collar job
t24 = 0.0756; % deviation of type 4 from type 1 in blue-collar job

% parameters for white-collar wage
% education effect
e11_true = 10.71;
e12_true = -0.17; % quadratic
% white-collar exp effect
e13_true = 10.45;
e14_true = -0.001; % quadratic
% blue-collar exp effect
e15_true = 2.18;


% parameters for blue-collar wage
% education effect
e21_true = 10.54;
e22_true = -0.22; % quadratic
% white-collar exp effect
e23_true = 2.29;
% blue-collar exp effect
e25_true = 14.7;
e26_true = -0.27; % quadratic

% parameters for school choice
e30_true = 4.89; % school intercept for type 1
t32 = 0.3352; % deviation of type 2 from type 1 in school choice
t33 = 0.0541; % deviation of type 3 from type 1 in school choice
t34 = -0.0226; % deviation of type 4 from type 1 in school choice

c1_true = 0.35; % tuition for college
c2_true = 3.92; % tuition for graduate school

% parameters for home choice
e40_true = 2.69;  % home intercept for type 1
t42 = -0.0215; % deviation of type 2 from type 1 in home choice
t43 = 1.6966; % deviation of type 3 from type 1 in home choice
t44 = 1.3128; % deviation of type 4 from type 1 in home choice

delta_true = 0.84; % discount rate

% parameters for var-cov of random shock
% we allow covariance between white- and blue-collar shocks
s11 = 0.153;
s22 = 0.34;
s33 = 3.82;
s44 = 4;
s12 = -0.04;

%}
% lump parameters into a list
param_true = [r1_true,r2_true, e11_true,e12_true, e13_true,e14_true,e15_true, ...
            e21_true,e22_true,e23_true,e25_true,e26_true, ...
            e30_true,c1_true,c2_true,e40_true,delta_true,s11,s22,s33,s44,s12, ...
            plow2,plow3,plow4,phigh2,phigh3,phigh4,...
            t12,t22,t32,t42,t13,t23,t33,t43,t14,t24,t34,t44];
param_est_m = [9.374177774	10.9031211	10.70406755	-0.170964288	10.52466794 ...
               -0.001004045	2.183049954	10.54349066	-0.220693643	2.300345331 ...
                14.72820152	-0.271384667	4.817715031	0.350739553	3.916424068 ...
                2.67666628	0.831698522	0.153481667	0.35004905	3.842698462 ...
           	    0.04	-0.040054487	0.364226901	0.501623437	0.090164157	...
                0.369223247	0.110141364	0.033089391	-0.599635302	0.140199499 ...
            	0.324625663	0.089994216	0.010017276	0.140255087	0.040054184 ...
            	0.200251729	0.390697579	0.280442679	-0.900934592	-0.500622793
];

param_est_f = [8.935277511 10.59691104 10.86953263 -0.173652606 9.376287833 ...
                -3.10E-05   2.15840407  10.42672116 -0.21070996  2.31459028 ...
                14.63985799  -0.278216076  4.549181363 0.350410702 3.956898787 ...
                4.170374769   0.833772691  0.154422703  0.344584113  3.83810966 ...
                4.9557        -4.00E-02  0.364226901	0.501623437	0.090164157	...
                0.369223247	0.110141364	0.033089391	-0.599635302	0.140199499 ...
            	0.324625663	0.089994216	0.010017276	0.140255087	0.040054184 ...
            	0.200251729	0.390697579	0.280442679	-0.900934592	-0.500622793];

[space,selectedspace] = select_space(fp, 10)
gooopy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           simulate initial education level          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random numbers for each element
rng("default")
randomNum = rand(fp.N, 1);

%%%%%%%%%%%% simulate initial education for males %%%%%%%%%%%
fp.ini_edu_male = zeros(fp.N, 1);

% Map random numbers to values based on probabilities
for i = 1:fp.N
    if randomNum(i) < fp.M7
        fp.ini_edu_male(i) = 7;
    elseif randomNum(i) < fp.M7+fp.M8
        fp.ini_edu_male(i) = 8;
    elseif randomNum(i) < fp.M7+fp.M8+fp.M9
        fp.ini_edu_male(i) = 9;
    elseif randomNum(i) < fp.M7+fp.M8+fp.M9+fp.M10
        fp.ini_edu_male(i) = 10;
    else
        fp.ini_edu_male(i) = 11;
    end
end

%%%%%%%% simulate initial education for females %%%%%%%%%%%
fp.ini_edu_female = zeros(fp.N, 1);

% Map random numbers to values based on probabilities
for i = 1:fp.N
    if randomNum(i) < fp.F7
        fp.ini_edu_female(i) = 7;
    elseif randomNum(i) < fp.F7+fp.F8
        fp.ini_edu_female(i) = 8;
    elseif randomNum(i) < fp.F7+fp.F8+fp.F9
        fp.ini_edu_female(i) = 9;
    elseif randomNum(i) < fp.F7+fp.F8+fp.F9+fp.F10
        fp.ini_edu_female(i) = 10;
    else
        fp.ini_edu_female(i) = 11;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          simulate data from true parameters         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
% create data
data = sim_funct(fp.N,fp,param_est_m);
%data = sim_funct(fp.N,fp,param_1);

% data moment
fp.data_mom = moment_funct_indirect(data,fp); % moments of the true parameters

fp.data_cell = cell;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  import real data                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlsy97male = readmatrix('E:\KW97new\KW1997 project\Backward Recursion Matlab\Formal Model\NLSY97male_new.csv');

nlsy97female = readmatrix('E:\KW97new\KW1997 project\Backward Recursion Matlab\Formal Model\NLSY97female_new.csv');


% data moment for male:
fp.male_mom = moment_funct_indirect(nlsy97male,fp); % moments of the true parameters

% data moment for female:
fp.female_mom = moment_funct_indirect(nlsy97female,fp); % moments of the true parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                estimate parameters                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% estimate nlsy97female
fp.option_est = optimset('LargeScale', 'off', 'Tolfun', 0.1, 'Tolx', 0.01, ...
    'Display', 'iter', 'MaxIter', 1e10, 'MaxFunEvals', 10000);
% estimation starting at true parameters
param_start = param_est_f;  

[param_est_female97, obj_eval_est, exitflag] = fminsearch(@(param)...
                obj_funct_indirect(param,fp,fp.female_mom), param_start, fp.option_est);

param_est_female97
goopy;
%}

%{
% estimate nlsy97male
fp.option_est = optimset('LargeScale', 'off', 'Tolfun', 0.001, 'Tolx', 0.01, ...
    'Display', 'iter', 'MaxIter', 1e10, 'MaxFunEvals', 10000);
% estimation starting at true parameters
param_start = param_est_m;  % note that param_true are the log of the true params

[param_est_male97, obj_eval_est, exitflag] = fminsearch(@(param)...
                obj_funct_indirect(param,fp,fp.male_mom), param_start, fp.option_est);

param_est_male97

goopy;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              test objective function                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
obj_eval = obj_funct_indirect(param_est_m,fp,fp.male_mom);

obj_eval

goopy;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        graph obj fct w.r.t one parameter            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this section is used for tuning parameters by hand

%{
obj_eval_list = [];
e11_grid = 7:0.1:12;
for j = 1:length(e11_grid)
param_start_here = param_est_m;
param_start_here(3) = e11_grid(j);
obj_eval_list(j) = obj_funct_indirect(param_start_here,fp,fp.male_mom)
end

[e11_grid' obj_eval_list'];
plot(e11_grid,obj_eval_list)
xlabel('e11');
goopy;



obj_eval_list = [];
r2_grid = 6:0.1:11;
for j = 1:length(r2_grid)
param_start_here = param_est_m;
param_start_here(2) = r2_grid(j);
obj_eval_list(j) = obj_funct_indirect(param_start_here,fp,fp.male_mom)
end

[r2_grid' obj_eval_list'];
plot(r2_grid,obj_eval_list)
xlabel('r2');
goopy;





obj_eval_list = [];
r1_grid = 6:0.1:11;
for j = 1:length(r1_grid)
param_start_here = param_est_m;
param_start_here(1) = r1_grid(j);
obj_eval_list(j) = obj_funct_indirect(param_start_here,fp,fp.male_mom)
end

[r1_grid' obj_eval_list'];
plot(r1_grid,obj_eval_list)
xlabel('r1');
goopy;





obj_eval_list = [];
delta_grid = 0.6:0.01:0.9;
for j = 1:length(delta_grid)
param_start_here = param_est_m;
param_start_here(17) = delta_grid(j);
obj_eval_list(j) = obj_funct_indirect(param_start_here,fp,fp.male_mom)
end

[delta_grid' obj_eval_list'];
plot(delta_grid,obj_eval_list)
xlabel('delta');
goopy;

obj_eval_list = [];
t44_grid = -2:0.1:2;
for j = 1:length(t44_grid)
param_start_here = param_1;
param_start_here(40) = t44_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t44_grid' obj_eval_list'];
plot(t44_grid,obj_eval_list)
xlabel('t44');
goopy;


obj_eval_list = [];
t34_grid = -2:0.1:2;
for j = 1:length(t34_grid)
param_start_here = param_1;
param_start_here(39) = t34_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t34_grid' obj_eval_list'];
plot(t34_grid,obj_eval_list)
xlabel('t34');
goopy;


obj_eval_list = [];
t24_grid = 0.6:0.01:2;
for j = 1:length(t24_grid)
param_start_here = param_1;
param_start_here(38) = t24_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t24_grid' obj_eval_list'];
plot(t24_grid,obj_eval_list)
xlabel('t24');
goopy;

obj_eval_list = [];
t14_grid = 0.6:0.1:3;
for j = 1:length(t14_grid)
param_start_here = param_1;
param_start_here(37) = t14_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t14_grid' obj_eval_list'];
plot(t14_grid,obj_eval_list)
xlabel('t14');
goopy;

obj_eval_list = [];
t43_grid = -2:0.1:2;
for j = 1:length(t43_grid)
param_start_here = param_1;
param_start_here(36) = t43_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t43_grid' obj_eval_list'];
plot(t43_grid,obj_eval_list)
xlabel('t43');
goopy;


obj_eval_list = [];
t33_grid = -0.2:0.01:0.3;
for j = 1:length(t33_grid)
param_start_here = param_1;
param_start_here(35) = t33_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t33_grid' obj_eval_list'];
plot(t33_grid,obj_eval_list)
xlabel('t33');
goopy;

obj_eval_list = [];
t23_grid = 0.8:0.01:1.2;
for j = 1:length(t23_grid)
param_start_here = param_1;
param_start_here(34) = t23_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t23_grid' obj_eval_list'];
plot(t23_grid,obj_eval_list)
xlabel('t23');
goopy;

obj_eval_list = [];
t13_grid = -0.2:0.01:0.2;
for j = 1:length(t13_grid)
param_start_here = param_1;
param_start_here(33) = t13_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t13_grid' obj_eval_list'];
plot(t13_grid,obj_eval_list)
xlabel('t13');
goopy;


obj_eval_list = [];
t42_grid = -2:0.05:0;
for j = 1:length(t42_grid)
param_start_here = param_1;
param_start_here(32) = t42_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t42_grid' obj_eval_list'];
plot(t42_grid,obj_eval_list)
xlabel('t42');

goopy;

obj_eval_list = [];
e40_grid = 0:0.05:6;
for j = 1:length(e40_grid)
param_start_here = param_1;
param_start_here(16) = e40_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[e40_grid' obj_eval_list'];
plot(e40_grid,obj_eval_list)
xlabel('e40');

goopy;



obj_eval_list = [];
t32_grid = 0.1:0.01:0.5;
for j = 1:length(t32_grid)
param_start_here = param_1;
param_start_here(31) = t32_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t32_grid' obj_eval_list'];
plot(t32_grid,obj_eval_list)
xlabel('t32');
goopy;

obj_eval_list = [];
t22_grid = 0:0.01:0.3;
for j = 1:length(t22_grid)
param_start_here = param_1;
param_start_here(30) = t22_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t22_grid' obj_eval_list'];
plot(t22_grid,obj_eval_list)
xlabel('t22');
goopy;

obj_eval_list = [];
t12_grid = -0.7:0.01:0.4;
for j = 1:length(t12_grid)
param_start_here = param_1;
param_start_here(29) = t12_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[t12_grid' obj_eval_list'];
plot(t12_grid,obj_eval_list)
xlabel('t12');
goopy;

obj_eval_list = [];
phigh4_grid = 0:0.01:0.5;
for j = 1:length(phigh4_grid)
param_start_here = param_1;
param_start_here(28) = phigh4_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[phigh4_grid' obj_eval_list'];
plot(phigh4_grid,obj_eval_list)
xlabel('phigh4');
goopy;

obj_eval_list = [];
phigh3_grid = 0:0.01:0.3;
for j = 1:length(phigh3_grid)
param_start_here = param_1;
param_start_here(27) = phigh3_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[phigh3_grid' obj_eval_list'];
plot(phigh3_grid,obj_eval_list)
xlabel('phigh3');
goopy;

obj_eval_list = [];
phigh2_grid = 0.2:0.01:0.5;
for j = 1:length(phigh2_grid)
param_start_here = param_1;
param_start_here(26) = phigh2_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[phigh2_grid' obj_eval_list'];
plot(phigh2_grid,obj_eval_list)
xlabel('phigh2');
goopy;

obj_eval_list = [];
plow4_grid = 0:0.01:0.4;
for j = 1:length(plow4_grid)
param_start_here = param_1;
param_start_here(25) = plow4_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[plow4_grid' obj_eval_list'];
plot(plow4_grid,obj_eval_list)
xlabel('plow4');
goopy;


obj_eval_list = [];
plow3_grid = 0:0.01:0.3;
for j = 1:length(plow3_grid)
param_start_here = param_1;
param_start_here(24) = plow3_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[plow3_grid' obj_eval_list'];
plot(plow3_grid,obj_eval_list)
xlabel('plow3');
goopy;

obj_eval_list = [];
plow2_grid = 0.2:0.01:0.5;
for j = 1:length(plow2_grid)
param_start_here = param_1;
param_start_here(23) = plow2_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.female_mom)
end

[plow2_grid' obj_eval_list'];
plot(plow2_grid,obj_eval_list)
xlabel('plow2');
goopy;


obj_eval_list = [];
t13_grid = -0.2:0.01:0.6;
for j = 1:length(t13_grid)
param_start_here = param_1;
param_start_here(33) = t13_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.male_mom)
end

[t13_grid' obj_eval_list'];
plot(t13_grid,obj_eval_list)
xlabel('t13');
goopy;


obj_eval_list = [];
t42_grid = -0.2:0.01:0.2;
for j = 1:length(t42_grid)
param_start_here = param_1;
param_start_here(32) = t42_grid(j);
obj_eval_list(j) = obj_funct(param_start_here,fp,fp.male_mom)
end

[t42_grid' obj_eval_list'];
plot(t42_grid,obj_eval_list)
xlabel('t42');
goopy;

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     bootstrap SE                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% predefine
boot_male_data_collect = nan(1018*fp.D, fp.data_cols, fp.B);

rand('state', 102);
draws = randi([1 1018], 1018, fp.B); % randomly draw 1018 id for bootstrap sample

fp.option_est = optimset('LargeScale', 'off', 'Tolfun', 0.1, 'Tolx', 0.01, ...
    'Display', 'iter', 'MaxIter', 1e10, 'MaxFunEvals', 10000);

param_est_boot = nan(fp.J,fp.B);

for b = 1:fp.B

    boot_data = [];

    for q = 1:1018
    % block of individuals with id == rand_order
    index = nlsy97male(:, fp.col_r) == draws(q, b);
    boot_data = [boot_data; nlsy97male(index,:)];
    end

    boot_male_data_collect(:,:,b) = boot_data;

    % data sample (created above)
    data_here = boot_male_data_collect(:,:,b);

    %moments for this data
    boot_mom = moment_funct(data_here,fp);
    
    %
    %estimate parameters using these moments
    param_start = param_est_m;
    
    [param_est_here,obj_eval_est,exitflag] = fminsearch(@(param) ...
        obj_funct(param,fp,boot_mom), ...
    param_start,fp.option_est);

    param_est_boot(:,b) = param_est_here';
    %}
end


boot_SE = std(param_est_boot')';

boot_SE

goopy;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         graph real data and simulated data          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
subplot(2,1,1),graph_funct(nlsy97male,fp,1018,fp.D)
xlabel('true data (male)');
subplot(2,1,2),graph_funct(data,fp,fp.N,fp.D)
xlabel('simulated data (male)');
goopy;
%}

%{
subplot(2,1,1),graph_funct(nlsy97female,fp,fp.N,fp.D)
xlabel('true data (female)');
subplot(2,1,2),graph_funct(data,fp,fp.N,fp.D)
xlabel('simulated data (female)');
goopy;
%}
















%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  bootstrap estimate data moment covariance matrix   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-define
boot_data_collect = nan(fp.N*fp.D, fp.data_cols, fp.B);
boot_moment = zeros(fp.M, 1, fp.E, fp.B);

% random number draws for bootstrap
% this creates fp.N x fp.B matrix of random data intergers
rand('state', 1131);
draws = randi([1 fp.N], fp.N, fp.B);

for b = 1:fp.B % fp.B = 200

    % randomly re-draw with replacement from original data
    boot_data = [];

    for q = 1:fp.N
    % block of individuals with id == rand_order
    index = nlsy97male(:, fp.col_r) == draws(q, b);
    boot_data = [boot_data; nlsy97male(index,:)];
    end
    
    boot_data_collect(:,:,b) = boot_data;
    
    % calculate moments for data draw
    boot_moment(:, b) = moment_funct(boot_data, fp);

end

% create bootstrapped cov matrix, off-diagonal elements set to zero

cov = zeros(fp.M, 1);
for b = 1:fp.B
    diff_here = boot_moment(:,b) - mean(boot_moment, 2);
    cov = cov + 1/fp.B.*diff_here.^2;
end

cov

%}
