close all; clear all; clc;
%%PART A
%%
% Question a
A = [ -0.0507   -3.861      0    -32.2
      -0.00117  -0.5164     1       0
      -0.000129  1.4168  -0.4932    0
          0        0        1       0];

B = [ 0
     -0.0717
     -1.645
      0];

% Calculate eigenvalues and eigenvectors
[V, D, W] = eig(A);

% Normalize left eigenvectors
W = W ./ diag(W' * V)';


eigenvalues = diag(D);
disp('Eigenvalues:');
disp(eigenvalues);


disp('Eigenvectors (right):');
disp(V);
disp('Eigenvectors (left):');
disp(W');

% Define initial conditions
x0_1 = [0; 0.1; 0; 0];
x0_2 = [0; -0.1; 0; 0];

% Compute the matrix exponential
expDt = diag(exp(eigenvalues));

% Compute x(t) for initial conditions
x_t_1 = V * expDt * inv(V) * x0_1;
x_t_2 = V * expDt * inv(V) * x0_2;

t_vec = linspace(0, 10, 1000);

x_t_1_num = zeros(4, length(t_vec));
x_t_2_num = zeros(4, length(t_vec));

% Compute the state trajectories 
for i = 1:length(t_vec)
    expDt_num = diag(exp(eigenvalues * t_vec(i)));
    x_t_1_num(:, i) = double(V * expDt_num * inv(V) * x0_1);
    x_t_2_num(:, i) = double(V * expDt_num * inv(V) * x0_2);
end

% Plot the results
state_names = {'change in speed', 'angle of attack', 'pitch rate', 'pitch angle'};

% Plot for x(0) = [0; 0.1; 0; 0]
fig1 = figure('Name', 'Initial Condition: x(0) = [0; 0.1; 0; 0]');
for i = 1:4
    subplot(2, 2, i);
    plot(t_vec, x_t_1_num(i, :));
    title([state_names{i}, ' for x(0) = [0; 0.1; 0; 0]']);
    xlabel('Time (s)');
    ylabel(state_names{i});
end

% Plot for x(0) = [0; -0.1; 0; 0]
fig2 = figure('Name', 'Initial Condition: x(0) = [0; -0.1; 0; 0]');
for i = 1:4
    subplot(2, 2, i);
    plot(t_vec, x_t_2_num(i, :));
    title([state_names{i}, ' for x(0) = [0; -0.1; 0; 0]']);
    xlabel('Time (s)');
    ylabel(state_names{i});
end

% X expressed as the modes of the system
exp_At = zeros(size(A));
for i = 1:length(eigenvalues)
    A_i = V(:, i) * W(:, i)'; 
    exp_At = exp_At + A_i * exp(eigenvalues(i));
end
x_t_1_modal =  exp_At * x0_1;
x_t_2_modal =  exp_At * x0_2;


%%
%%

% Question b
T = 4; % sampling period
t = 0:0.01:10; 

u_b = -1*(t>= 0 & t <= T);
sys = ss(A, B, eye(4), 0);  

x0 = [0; 0; 0; 0];
x = lsim(sys, u_b, t, x0'); 

fig3 = figure('Name', 'Initial Condition: x(0) = [0; 0; 0; 0] and u = -1 for 0=<t<=T');
for i = 1:4
    subplot(2, 2, i);
    plot(t, x(:, i));
    title(state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end


%% question c)
u_c = -1 * (t>=0);  % negative unit step input (u = -1 for t >= 0)
x_c = lsim(sys, u_c, t, x0); % again y(t)=x(t)

% plot with subplots
fig4 = figure('Name', 'Initial Condition: x(0) = [0; 0; 0; 0] and u = -1 for 0=<t<=T');
for i = 1:4
    subplot(2, 2, i);
    plot(t, x_c(:, i));
    title(state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end

C = [0 0 1 0]; %Cause of question z (y=q=x3=pitch rate)
%Question e
rank(ctrb(A,B))
%Question z
rank(obsv(A,C))








































