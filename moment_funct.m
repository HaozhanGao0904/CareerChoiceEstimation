function mom = moment_funct(data, fp)

global cell

t_vec = data(:, fp.col_t);
d_vec = data(:, fp.col_d);
x1_vec = data(:, fp.col_x1);
x2_vec = data(:, fp.col_x2);
g_vec = data(:, fp.col_g);

mom = nan(fp.M, 1, fp.E);
cell = mom;

% choice fractions conditional on all possible states

for k = 1:fp.E % loop over initial education levels
    i = 1;
    ini_g = 6+k; % initial education level
% there are 1002 states and 12024 moments
    for t = 11:-1:1 % for t= 10,9,8,...,1
    
        for x1_index = 1:t % loop over white-collar work exp states
            x1 = x1_index - 1; % actual white-collar work exp
            
            for x2_index = 1: t-x1 % loop over blue-collar work exp states
                x2 = x2_index - 1; % actual blue-collar work exp
    
                for g_index = 1: t-x1-x2 % loop over education years states
                    g = (g_index - 1) + ini_g; % actual education years
    
                    % for each state, store the folloing information:
                    j = 12*(i-1);
                    
                    for d = 1: 4 % for each choice in each state
    
                        totalnum = nansum(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g);
    
                        % store the choice fractions in moments
                        mom(j+(d-1)*3+1,1,k) = nansum(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g & d_vec == d)/totalnum;
                        cell(j+(d-1)*3+1,1,k) = nansum(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g & d_vec == d);
    
                        % store the avg. wage in moments
                        mom(j+(d-1)*3+2,1,k) = nanmean(data(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g, fp.col_lnw))/100;
                        cell(j+(d-1)*3+2,1,k) = nansum(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g & d_vec == d)/10;
                        
                        % store the std wages in moments
                        mom(j+(d-1)*3+3,1,k) = nanstd(data(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g, fp.col_lnw))/100;
                        cell(j+(d-1)*3+3,1,k) = nansum(t_vec == t & x1_vec == x1 & x2_vec == x2 & g_vec == g & d_vec == d)/10;
    
                    end
                                
                    i = i+1;
                
                end
            end
        end
    end


    for i = 1:fp.M
        if isnan(mom(i,1,k))
            mom(i,1,k) = 0;
        end
    end
end










