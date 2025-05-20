clear all; close all; clc;
%%
% Question α
A = [ -0.0507   -3.861      0    -32.2
      -0.00117  -0.5164     1       0
      -0.000129  1.4168  -0.4932    0
          0        0        1       0];

B = [ 0
     -0.0717
     -1.645
      0];
state_names = {'change in speed', 'angle of attack', 'pitch rate', 'pitch angle'};
controller_eigenvalues = [-1.25 + 2.2651j, -1.25 - 2.2651j, -0.01 + 0.095j, -0.01 - 0.095j];

K = place(A,B,controller_eigenvalues);

disp('Feedback gain matrix K:');
disp(K);

%%
%Question β
x0_1 = [0; 0.1; 0; 0];
x0_2 = [0; -0.1; 0; 0];

%new system
% x' = Ax= Bu = Ax -BKx = (A-BK)x

T = 4; % sampling period
t = 0:0.01:400; 

r = zeros(size(t)); %arbitrary value
sys = ss(A-B*K, B, eye(4), 0);  

x_1 = lsim(sys,r,t,x0_1);
x_2 = lsim(sys,r,t,x0_2);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, x_1(:, i));
    title("Closed loop system for "+ "x(0)=[0, 0.1, 0, 0]" +":",state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end


figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, x_2(:, i));
    title("Closed loop system for "+ "x(0)=[0, -0.1, 0, 0]" +":",state_names{i});    xlabel('Time (s)');
    ylabel(state_names{i});
end


%%
%Question γ
x0 = [0; 0; 0; 0];
r_b_1 = -1*(t>= 0 & t <= T);
xg_1 = lsim(sys,r_b_1,t,x0_1);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, xg_1(:, i));
    title("Closed loop system for "+ "x(0)=[0, 0, 0, 0]" + "and r= -1 for 0=<t<=T" +":",state_names{i});    xlabel('Time (s)');
    xlabel('Time (s)');
    ylabel(state_names{i});
end

r_b_2 = -1 * (t>=0);
xg_2 = lsim(sys,r_b_2,t,x0_1);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, xg_2(:, i));
    title("Closed loop system for "+ "x(0)=[0, 0, 0, 0]" + "and r= -1" +":",state_names{i});    xlabel('Time (s)');
    xlabel('Time (s)');
    ylabel(state_names{i});
end
%%
%Question δ

C = [0 0 1 0];
observer_eigenvalues = [-0.1, -0.421, -0.587, -1];
L = place(A', C', observer_eigenvalues)';

disp('Observer gain matrix L:');
disp(L);

%%
%Question ε
%e(t) = x(t) - x_hat(t)
%e(t)' = (A-LC)e(t)
%e(t) = exp(A-LC)*e(0)

x0_hat = [0.2; -0.1; 0.1; -0.1];

e0_1 = x0_1 - x0_hat;
e0_2 = x0_2 - x0_hat;

sys_e = ss(A-L*C,B,eye(4),0);

e1 = lsim(sys_e,r,t,e0_1);
e2 = lsim(sys_e,r,t,e0_2);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, e1(:, i));
    title("Error of:",state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, e2(:, i));
    title("Error of:",state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end

%%
%Question ζ
% x_bar = (x,e)^T
t_z = 0:0.01:400; 
r_z = zeros(size(t_z));

A_bar = [A - B * K, B * K; zeros(4, 4), A - L * C];

x0_bar1 = [x0_1; e0_1];
x0_bar2 = [x0_2; e0_2];

sys_bar = ss(A_bar, zeros(8, 1), eye(8), zeros(8, 1));

x_bar1 = lsim(sys_bar, r_z, t_z, x0_bar1);
x_bar2 = lsim(sys_bar, r_z, t_z, x0_bar2);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t_z, x_bar1(:, i),"b");
    hold on;
    plot(t_z, -x_bar1(:, i+4)+x_bar1(:, i),"r");
    title([state_names{i}, ' and ', "estimation"]);
    xlabel('Time (s)');
    ylabel(state_names{i});
end

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t_z,x_bar2(:, i),"b");
    hold on;
    plot(t_z, -x_bar2(:, i+4)+x_bar2(:, i),"r");
    title([state_names{i}, ' and ', "estimation"]);
    xlabel('Time (s)');
    ylabel('State and Error');
    legend("Actual Value", "Estimation")
end

%%
%Question η

B_bar = [ B
       zeros(size(B))];

C_bar = [ C    zeros(size(C))];

D = zeros(8, 1);
new_sys = ss(A_bar, B_bar, C_bar, 0);

systf = tf(new_sys);

rank(ctrb(A_bar, B_bar))

e_eig = eig(A-L*C);

disp("Eigenvalues of e");
disp(e_eig);

%%
%Question θ

Q = [ 10  0  0  0       
      0  10 0  0
      0  0  100 0
      0  0  0  100];

R=1;

K_lqr = lqr(A,B,Q,R);

disp("Optimal K");
disp(K_lqr);

sys_lqr = ss(A-B*K_lqr,B,eye(4),0);


[V_lqr,eig_lqp, W_lqr] = eig(A-B*K_lqr);

disp("Eigvalues of system");
eig_lqp1 = diag(eig_lqp)
disp(eig_lqp1);


x_lqr1 = lsim(sys_lqr, r_z, t_z, x0_1);
x_lqr2 = lsim(sys_lqr, r_z, t_z, x0_2);

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, x_lqr1(:, i));
    title("Closed loop system with optimal K for "+ "x(0)=[0, 0.1, 0, 0]" +":",state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(t, x_lqr2(:, i));
    title("Closed loop system with optimal K for "+ "x(0)=[0, 0.1, 0, 0]" +":",state_names{i});
    xlabel('Time (s)');
    ylabel(state_names{i});
end


sys_K = ss(A-B*K,B,C,0);
sys_Klqr = ss(A-B*K_lqr,B,C,0);

x1_K = lsim(sys_K, r_z, t_z, x0_1);
x2_K = lsim(sys_K, r_z, t_z, x0_2);

x1_Klqr = lsim(sys_Klqr, r_z, t_z, x0_1);
x2_Klqr = lsim(sys_Klqr, r_z, t_z, x0_2);

figure;
subplot(2, 1, 1);
plot(t, x1_K);
title("x3 with K with pole placement "+ "x(0)=[0, 0.1, 0, 0]")
xlabel('Time (s)');

subplot(2, 1, 2);
plot(t, x1_Klqr);
title("x3 with K with LQR "+ "x(0)=[0, 0.1, 0, 0]")
xlabel('Time (s)');

figure;
subplot(2, 1, 1);
plot(t, x2_K);
title("x3 with K with pole placement "+ "x(0)=[0, -0.1, 0, 0]")
xlabel('Time (s)');

subplot(2, 1, 2);
plot(t, x2_Klqr);
title("x3 with K with LQR "+ "x(0)=[0, -0.1, 0, 0]")
xlabel('Time (s)');





