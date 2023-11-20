function data_space = extract_space(gender)

if gender == "m"
data = readmatrix('Z:\KW1997 project\Backward Recursion Matlab\Formal Model\NLSY97male_lnwage.csv');
else
data = readmatrix('Z:\KW1997 project\Backward Recursion Matlab\Formal Model\NLSY97female_lnwage.csv');
end

data_space = {};
data_space_t = {};
t_vec = data(:, 2);
x1_vec = data(:, 4);
x2_vec = data(:, 5);
g_vec = data(:, 6);


for t = 1:11
indices_t = (t_vec == t); % create a boolean array to extract time = t

% extract x1, x2, g vectors when time = t
x1_t = x1_vec(indices_t);
x2_t = x2_vec(indices_t);
g_t = g_vec(indices_t);

unique_cells = {};  % Initialize a cell array to store unique cells

for i = 1:length(x1_t)
    current_cell = [x1_t(i), x2_t(i), g_t(i)];

    % Check if the current cell already exists in unique_cells
    is_duplicate = false;
    for j = 1:length(unique_cells)
        if isequal(current_cell, unique_cells{j})
            is_duplicate = true;
            break;
        end
    end

    if ~is_duplicate
        data_space_t{end+1} = current_cell;
        unique_cells{end+1} = current_cell;
    end
end

data_space{t} = data_space_t;
end




