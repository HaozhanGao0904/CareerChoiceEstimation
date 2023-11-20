% objective function

function obj_eval = obj_funct(param,fp,data_mom)
global g_diff

% 0) check if the cov matrix for random shocks are positive semi-definite
% if not, kill the function with a large objective function value

s11 = param(18);
s22 = param(19);
s12 = param(22);

cov_matrix = [s11, s12; s12, s22;];
evalues = eig(cov_matrix);

if any(evalues < 0)
    obj_eval = 1e10;
    return;
end

% 1) simulate data
sim_data = sim_funct(fp.N,fp,param);


% 2) calculate moments from simulated data
sim_mom = moment_funct(sim_data,fp);

% 3) calculate g_diff
g_diff = nan(fp.M,1,fp.E); % g_diff is the difference between the simulated moments and the real moments

for k = 1:fp.E
    g_diff(:,1,k) = data_mom(:,1,k) - sim_mom(:,1,k);
end

% 4) form objective function

obj_eval = sum(g_diff(:,1,1).*g_diff(:,1,1)) + sum(g_diff(:,1,2).*g_diff(:,1,2)) + sum(g_diff(:,1,3).*g_diff(:,1,3)) + ...
           sum(g_diff(:,1,4).*g_diff(:,1,4)) + sum(g_diff(:,1,5).*g_diff(:,1,5));