
% this function will create a list of state spaces for every time period
function [space,selectedspace] = select_space(fp, ini_edu_index)
fp.T = 50;
interp_num = fp.Interp;
space = {};
selectedspace = {};

% first we store all states in "space"

for t = fp.T:-1:2 % loop over time period
    spaceofT = {};
    count = 1;
    for x1_index = 1:t
        x1 = x1_index - 1;
        for x2_index = 1: t-x1 % loop over blue-collar work exp
            x2 = x2_index - 1;
            for g_index = 1: t-x1-x2 % loop over education
            g = 6 + ini_edu_index + (g_index - 1); % actual education years
            spaceofT{count} = [x1,x2,g];
            count = count + 1;
            end
        end
    end
    space{t} = spaceofT;
end

% then we select a subset of space that is used for interpolation


for t = fp.T:-1:2
    len = length(space{t});
    if len > interp_num
        rng("default");
        indices = randperm(len, interp_num);
        selectedspace{t} = space{t}(indices(1:interp_num));
    else
        selectedspace{t} = space{t};
    end
   
end



