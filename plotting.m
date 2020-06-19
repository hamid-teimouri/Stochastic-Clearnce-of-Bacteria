% ===============================================================================
% Testing Simulation of Stochastic Clearance of Bacteria with an exact solution for N=2
% Dr. Hamid Teimouri, Rice University - Department of Chemistry  
% ===============================================================================
close all; clc;

lambda1 = 3.0/60;
lambda2 = 0.3/60; 
lambda=(lambda1 + lambda2)/2;
gama = lambda;
delta = gama;
a = (lambda1 - lambda2)/2;

%% upload data from simulation
Array1=textread('f_MC.txt');
Array2=textread('T_MC.txt');
Array3=textread('f1_MC.txt');
Array4=textread('f2_MC.txt');
Array5=textread('T1_MC.txt');
Array6=textread('T2_MC.txt');

x1 = Array1(:,1); f_mc = Array1(:,2); 

x2 = Array2(:,1); T_mc = Array2(:,2);

x3 = Array3(:,1); f1_mc = Array3(:,2);

x4 = Array4(:,1); f2_mc = Array4(:,2);

x5 = Array5(:,1); T1_mc = Array5(:,2);

x6 = Array6(:,1); T2_mc = Array6(:,2);

%% define vectors for Exact results
f1_exact = zeros(length(x1),1);
f2_exact = zeros(length(x1),1);
T1_exact = zeros(length(x1),1);
T2_exact = zeros(length(x1),1);

% Exact solution of the moded for N=2
for i = 1:length(x1)
    phi(i) = lambda2* x1(i);
    x(i)=x1(i);
    
    % Exact extinction probabilities:
    f1_exact(i) = (phi(i) * (2*delta + lambda2 + phi(i)))/((gama + lambda1 + phi(i))*(delta + lambda2 + phi(i)) - gama*delta);
    f2_exact(i) = (phi(i) * (2.0*gama + lambda1 + phi(i)))/((gama + lambda1 + phi(i))*(delta + lambda2 + phi(i)) - gama*delta);
    f_exact(i)=x(i)./(1 + x(i));
    
    % Exact probability ratio r_f:
    rf_exact(i) = (f1_exact(i)  + f2_exact(i) )./(2*f_exact(i));
   
    %========================================================
    % Exact extinction times:
    T1_num(i) = 2*gama*delta + 2*(delta^2) + delta*(lambda1 + 3*lambda2 + 4*phi(i)) + (lambda2 + phi(i))^2;
    T1_den(i) = (2*delta + lambda2 + phi(i))*( (lambda2 + phi(i))*(gama + lambda1 + phi(i)) + delta*(lambda1 + phi(i)) );
    T1_exact(i) = T1_num(i)./T1_den(i);

    T2_num(i) = 2*(gama^2) + gama*(2*delta + 3*lambda1 + lambda2 + 4*phi(i)) + (lambda1 + phi(i))^2;
    T2_den(i) = (2*gama + lambda1 + phi(i))*( (lambda2 + phi(i))*(gama + lambda1 + phi(i)) + delta*(lambda1 + phi(i)) );
    T2_exact(i) = T2_num(i)./T2_den(i);
    T_exact(i)=1/((1 + x(i))*lambda);

    % Exact time ratio r_T     
    rT_exact(i) = (T1_exact(i)  + T2_exact(i) )./(2*T_exact(i));
     
end

rf_mc = (f1_mc + f2_mc)./(2*f_mc);
rT_mc = (T1_mc + T2_mc)./(2*T_mc);


%% Plots
% test results for extiction probability
figure()
subplot(2,2,1)
plot(x3,f1_mc,'-',x3,f1_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('f_1', 'FontSize', 20);

subplot(2,2,2)
plot(x3,f2_mc,'-',x3,f2_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('f_2', 'FontSize', 20);

subplot(2,2,3)
plot(x3,f_mc,'-',x3,f_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('f', 'FontSize', 20);

subplot(2,2,4)
plot(x1,rf_mc,'-',x1,rf_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('r_f', 'FontSize', 20);

sgtitle('Comparison of exact and MC simulation results for extinction probability')

%==========================================
% test results for extiction times
figure()
subplot(2,2,1)
plot(x3,T1_mc,'-',x3,T1_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('T_1', 'FontSize', 20);

subplot(2,2,2)
plot(x3,T2_mc,'-',x3,T2_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('T_2', 'FontSize', 20);

subplot(2,2,3)
plot(x3,T_mc,'-',x3,T_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('T', 'FontSize', 20);

subplot(2,2,4)
plot(x3,rT_mc,'-',x3,rT_exact,'d','LineWidth',2)
legend('Simulation','Exact solution')
xlabel('x', 'FontSize', 15);
ylabel('r_T', 'FontSize', 20);

sgtitle('Comparison of exact and MC simulation results for extinction time')
