function graph = graph_funct(data,fp,popnum, t)

global cell

t_vec = data(:, fp.col_t);
d_vec = data(:, fp.col_d);
x1_vec = data(:, fp.col_x1);
x2_vec = data(:, fp.col_x2);
g_vec = data(:, fp.col_g);

d1 = zeros(1,t);
d2 = zeros(1,t);
d3 = zeros(1,t);
d4 = zeros(1,t);

for i = 1:t
    d1(i) = nansum(t_vec == i & d_vec == 1)/popnum;
    d2(i) = nansum(t_vec == i & d_vec == 2)/popnum;
    d3(i) = nansum(t_vec == i & d_vec == 3)/popnum;
    d4(i) = nansum(t_vec == i & d_vec == 4)/popnum;
end

% Plot all vectors together
time = 1:t;
plot(time, d1, 'LineWidth', 2);
hold on;  % This allows adding more plots to the same figure
plot(time, d2, 'LineWidth', 2);
plot(time, d3, 'LineWidth', 2);
plot(time, d4, 'LineWidth', 2);
hold off; % Stop adding to the same figure
grid on;
%title('career choices over time');
xlabel('Time Period');
ylabel('choice fraction');
legend('white-collar', 'blue-collar', 'school', 'home');