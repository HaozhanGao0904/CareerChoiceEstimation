function sim_mom = moment_funct_indirect(data, fp)

    warning('off', 'MATLAB:rankDeficientMatrix'); % suppress the warning for regressions in the moment function

    % Extract data columns based on provided indices
    t_vec = data(:, fp.col_t);
    d_vec = data(:, fp.col_d);
    x1_vec = data(:, fp.col_x1);
    x2_vec = data(:, fp.col_x2);
    g_vec = data(:, fp.col_g);
    lnw_vec = data(:, fp.col_lnw);

    % Initialize the coefficients matrix
    sim_mom_matrix = NaN(20,11); 

    % Loop over each period
    for t = 1:11
        % Extract data for the current period
        current_indices = t_vec == t;
        
        % Probit regressions for each choice
        for choice = 1:4
            dep_var = d_vec(current_indices) == choice;
            
            try
                probit_result = glmfit(indep_vars, dep_var, 'binomial', 'link', 'probit');
                sim_mom_matrix((choice-1)*4+1:(choice-1)*4+4, t) = probit_result;
            catch
                sim_mom_matrix((choice-1)*4+1:(choice-1)*4+4, t) = nan;
            end
        
        end

       % OLS regression for wage
        wage_indices = current_indices & ~isnan(lnw_vec);
        if any(wage_indices)
            dep_var_lnw = lnw_vec(wage_indices);
            indep_vars_lnw = [ones(sum(wage_indices), 1), g_vec(wage_indices), x1_vec(wage_indices), x2_vec(wage_indices)];
            try
                ols_result = regress(dep_var_lnw, indep_vars_lnw);
                sim_mom_matrix(17:20, t) = ols_result;
            catch
                sim_mom_matrix(17:20, t) = nan;
            end
        else
            sim_mom_matrix(17:20, t) = nan;
        end
    end

    % Flatten the matrix into a vector
    sim_mom = sim_mom_matrix(:);
end
