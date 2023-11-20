
% this function will create a list of state spaces for every time period
%function space = interpolation_test(fp, ini_edu)
T = 50;
interp_num = 200;
space = {};
selectedspace = {};

% first we store all states in "space"

for t = T:-1:2 % loop over time period
    spaceofT = {};
    count = 1;
    for x1_index = 1:t
        x1 = x1_index - 1;
        for x2_index = 1: t-x1 % loop over blue-collar work exp
            x2 = x2_index - 1;
            for g_index = 1: t-x1-x2 % loop over education
            g = 6 + (g_index - 1); % actual education years
            spaceofT{count} = [x1,x2,g];
            count = count + 1;
            end
        end
    end
    space{t} = spaceofT;
end

% then we select a subset of space that is used for interpolation
for t = T:-1:2
    selectedspaceoft = {};
    len = length(space{t});

    if len <= interp_num
        selectedspaceoft = space{t};
    else
        stepsize = floor(len/interp_num);
        i = 1;
        j = 1;
        while i <= len
            selectedspaceoft{j} = space{t}{i};
            i = i + stepsize;
            j = j + 1;
        end
    end
    selectedspace{t} = selectedspaceoft;
end

% we simulate EMAX for selected space for each t and run regression to get
% the coeffiencts

for t = T:-1:2
    get_EMAX
end

