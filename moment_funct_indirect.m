function sim_mom = moment_funct_indirect(data, fp)

    warning('off', 'MATLAB:rankDeficientMatrix'); % suppress the warning for regressions in the moment function

    % Extract data columns based on provided indices
    t_vec = data(:, fp.col_t);
    d_vec = data(:, fp.col_d);
    x1_vec = data(:, fp.col_x1);
    x2_vec = data(:, fp.col_x2);
    g_vec = data(:, fp.col_g);
    lnw_vec = data(:, fp.col_lnw);

    % Initialize the coefficients matrix with 30 rows for the probit and OLS regressions
    sim_mom_matrix = NaN(30, 11); % Adjusted size to 30 x 11 for all coefficients

    % Loop over each period
    for t = 1:11
        % Extract data for the current period
        current_indices = t_vec == t;

        % Probit regressions for each choice, including x1^2 and x2^2
        for choice = 1:4
            dep_var = d_vec(current_indices) == choice;
            % Include the squares of x1 and x2 as additional regressors
            indep_vars = [ones(sum(current_indices), 1), ...
                          g_vec(current_indices), ...
                          x1_vec(current_indices), ...
                          x2_vec(current_indices), ...
                          x1_vec(current_indices).^2, ...
                          x2_vec(current_indices).^2];
            try
                % glmfit returns a vector of coefficients
                probit_result = glmfit(indep_vars, dep_var, 'binomial', 'link', 'probit');
                % Store the probit result coefficients
                sim_mom_matrix((choice-1)*6+1:choice*6, t) = probit_result;
            catch
                sim_mom_matrix((choice-1)*6+1:choice*6, t) = nan;
            end
        end

        % OLS regression for wage, including x1^2 and x2^2
        wage_indices = current_indices & ~isnan(lnw_vec);
        if any(wage_indices)
            dep_var_lnw = lnw_vec(wage_indices);
            % Include the squares of x1 and x2 as additional regressors
            indep_vars_lnw = [ones(sum(wage_indices), 1), ...
                              g_vec(wage_indices), ...
                              x1_vec(wage_indices), ...
                              x2_vec(wage_indices), ...
                              x1_vec(wage_indices).^2, ...
                              x2_vec(wage_indices).^2];
            try
                % regress returns a vector of coefficients
                ols_result = regress(dep_var_lnw, indep_vars_lnw);
                % Store the OLS result coefficients
                sim_mom_matrix(25:30, t) = ols_result;
            catch
                sim_mom_matrix(25:30, t) = nan;
            end
        else
            sim_mom_matrix(25:30, t) = nan;
        end
    end

    % Flatten the matrix into a vector
    sim_mom = sim_mom_matrix(:);
end
