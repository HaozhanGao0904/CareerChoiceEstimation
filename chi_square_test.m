function [chi2stats, critical_value] = chi_square_test(real_data, simulated_data)
    num_periods = 11; % Number of periods
    num_choices = 4; % Number of occupation choices
    chi2stats = zeros(num_periods, num_choices); % Store Chi-square stats for each period
    critical_value = chi2inv(0.95, num_choices - 1); % Common critical value for all periods

    for period = 1:num_periods
        % Calculate observed frequencies for the real data
        observed = zeros(1, num_choices);
        for choice = 1:num_choices
            observed(choice) = sum(simulated_data(:,2) == period) * sum(real_data(:,2) == period & real_data(:,3) == choice)...
                            /sum(real_data(:,2) == period);
        end
        
        % Calculate expected frequencies for the simulated data
        expected = zeros(1, num_choices);
        for choice = 1:num_choices
            expected(choice) = sum(simulated_data(:,2) == period & simulated_data(:,3) == choice);
        end
        
        % Ensure no zero observed frequencies to avoid division by zero
        observed(observed == 0) = eps;
        
        % Calculate the Chi-square statistic for each choice
        chi2stat = ((observed - expected).^2) ./ observed;

        % Store the Chi-square statistics for this period
        chi2stats(period, :) = chi2stat;
    end
end
