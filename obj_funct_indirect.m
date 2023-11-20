function obj_eval = your_function_name(param, fp, data_mom)

    % 0) Check if the covariance matrix for random shocks is positive semi-definite
    s11 = param(18);
    s22 = param(19);
    s12 = param(22);

    cov_matrix = [s11, s12; s12, s22];
    evalues = eig(cov_matrix);

    if any(evalues < 0)
        obj_eval = 1e10;
        return;
    end

    % 1) Simulate data
    sim_data = sim_funct(fp.N, fp, param);

    % 2) Calculate moments from simulated data
    sim_mom = moment_funct_indirect(sim_data, fp);

    % 3) Calculate g_diff
    g_diff_matrix = data_mom - sim_mom;
    g_diff_matrix(isnan(g_diff_matrix)) = 0; % convert nans to 0s
    obj_eval = g_diff_matrix' * g_diff_matrix; % This computes the squared differences in a matrix form
    obj_eval = 100*sum(obj_eval(:)); % This sums up all the elements of the matrix to get a scalar value

end
