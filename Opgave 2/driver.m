%% Exercise 1
clc;
clear all;
fs = 13;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
%% Task 2
n = 10; u = 0.2; d = 1;
[H, g, A, b] = ConstructEqQP(n,u,d,1);

% Task3
[K,h] = KKTSystem(n,u,d,1);

%% Task 4
n_used=length(10:25:1000);
t_elapse = zeros(4,n_used);
avg = 40;

i=0;
for n=10:25:1000
    i=i+1;
    
    [K,h] = KKTSystem(n,u,d,0); % choose sparsity
    [H, g, A, b] = ConstructEqQP(n,u,d,0); % choose sparsity
    
    tic;
    for ii=1:avg
        [~,~] = KKTLUSolve(n,H,g,A,b,K,h);
    end
    e=toc/avg;
    t_elapse(1,i) = e;
    
    tic;
    for ii=1:avg
        [~,~] = KKTLDLSolve(n,H,g,A,b,K,h);
    end
    e1 = toc/avg;
    t_elapse(2,i) = e1;
    
    tic;
    for ii=1:avg
        [~,~] = KKTNSSolve(n,H,g,A,b,K,h);
    end
    e2 = toc/avg;
    t_elapse(3,i) = e2;
    
    tic;
    for ii=1:avg
        [~,~] = KKTRSSolve(n,H,g,A,b,K,h);
    end
    e3 = toc/avg;
    t_elapse(4,i) = e3;
end

%% Sparse
t_elapse2 = zeros(4,n_used);
avg = 80;

i=0;
for n=10:25:1000
    i=i+1;
    
    [K,h] = KKTSystem(n,u,d,1); % choose sparsity
    [H, g, A, b] = ConstructEqQP(n,u,d,1); % choose sparsity
    
    tic;
    for ii=1:avg
        [~,~] = KKTLUSolve(n,H,g,A,b,K,h);
    end
    e=toc/avg;
    t_elapse2(1,i) = e;
    
    tic;
    for ii=1:avg
        [~,~] = KKTLDLSolve(n,H,g,A,b,K,h);
    end
    e1 = toc/avg;
    t_elapse2(2,i) = e1;
end

%% CPU Time
t_elapse3 = zeros(4,n_used);
avg = 40;
i=0;
for n=10:25:1000
    i=i+1;
    [K,h] = KKTSystem(n,u,d,0); % choose sparsity
    [H, g, A, b] = ConstructEqQP(n,u,d,0); % choose sparsity
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTLUSolve(n,H,g,A,b,K,h);
    end
    e=(cputime-t)/avg;
    t_elapse3(1,i) = e;
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTLDLSolve(n,H,g,A,b,K,h);
    end
    e1 = (cputime-t)/avg;
    t_elapse3(2,i) = e1;
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTNSSolve(n,H,g,A,b,K,h);
    end
    e2 = (cputime-t)/avg;
    t_elapse3(3,i) = e2;
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTRSSolve(n,H,g,A,b,K,h);
    end
    e3 = (cputime-t)/avg;
    t_elapse3(4,i) = e3;
end

%%
t_elapse4 = zeros(2,n_used);
avg = 800;
i=0;
for n=10:25:1000
    i=i+1;
    
    [K,h] = KKTSystem(n,u,d,1); % choose sparsity
    [H, g, A, b] = ConstructEqQP(n,u,d,1); % choose sparsity
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTLUSolve(n,H,g,A,b,K,h);
    end
    e = (cputime-t)/avg;
    t_elapse4(1,i) = e;
    
    t = cputime;
    for ii=1:avg
        [~,~] = KKTLDLSolve(n,H,g,A,b,K,h);
    end
    e1 = (cputime-t)/avg;
    t_elapse4(2,i) = e1;
end

%%

figure(1);
n1 = 10:25:1000;
ax1=subplot(1,2,1);
hold on;
set(gca, 'YScale', 'log')
plot(n1,t_elapse(1,:),'Linewidth',2)
plot(n1,t_elapse(2,:),'Linewidth',2)
plot(n1,t_elapse(3,:),'Linewidth',2)
plot(n1,t_elapse(4,:),'Linewidth',2)
plot(n1,t_elapse2(1,:),'Linewidth',2)
plot(n1,t_elapse2(2,:),'Linewidth',2)
grid on
hold off;
set(gca,'FontSize',14)
%title('Performance of QP-solver vs. problem size','Interpreter','Latex','Fontsize',17);
legend(ax1,{'LU','LDL','Null-space', 'Range-space','LU-Sparse','LDL-Sparse'},'Interpreter','Latex','Fontsize',12, 'Location', 'northwest');
xlabel(ax1,'$n$','Interpreter','Latex','Fontsize',17)
ylabel(ax1,'Average wall-time','Interpreter','Latex','Fontsize',17)

ax2=subplot(1,2,2);
hold on;
set(gca, 'YScale', 'log')
plot(n1,t_elapse3(1,:),'Linewidth',2)
plot(n1,t_elapse3(2,:),'Linewidth',2)
plot(n1,t_elapse3(3,:),'Linewidth',2)
plot(n1,t_elapse3(4,:),'Linewidth',2)
plot(n1,t_elapse4(1,:),'Linewidth',2)
plot(n1,t_elapse4(2,:),'Linewidth',2)
grid on
hold off;
set(gca,'FontSize',14)
title('Performance of QP-solver vs. problem size','Interpreter','Latex','Fontsize',17);
legend({'LU','LDL','Null-space', 'Range-space','LU-Sparse','LDL-Sparse'},'Interpreter','Latex','Fontsize',12, 'Location', 'northwest');
xlabel('$n$','Interpreter','Latex','Fontsize',17)
ylabel('Average cpu-time','Interpreter','Latex','Fontsize',17)
set(get(gca,'title'),'Position',[-220 1.2 15])

%% Sparsity pattern
n=100
[H, g, A, b] = ConstructEqQP(n,u,d,1);

% Task3
[K,h] = KKTSystem(n,u,d,1);

spy(K, 8)
title('Sparsity of KKT matrix','Interpreter','Latex','Fontsize',17)
xlabel('')
%% Exercise 2: 
A = [1,0,1,1,-5; 0, 1 ,-1,-5,1];
b= [0,0,-2,-20,-15]';
g = [1 , -2];
f=@(x1,x2) (x1-1).^2 + (x2-2.5).^2;


x = -5:0.05:5;
y = -5:0.05:5;
[X,Y] = meshgrid(x,y);
f=@(x1,x2) (x1-1).^2 + (x2-2.5).^2;
F = f(X,Y);
[c,h] = contour(X,Y,F,100,'linewidth',3);
colorbar

yc1 = 1+1/2*x;
yc2 = 3-1/2*x;
yc3= 1/2*x-1;
yc4 = zeros(length(x));


hold on
fill(x,yc1, [0.2, 0.2, 0.2])
fill(x,yc2, [0.2, 0.2, 0.2])
fill(x,yc3, [0.2, 0.2, 0.2])
fill(x, yc4, [0.2, 0.2, 0.2])
fill(yc4, x, 'r')
hold off;
xlim([-5,5])
ylim([-5,5])
xlabel('x1')
ylabel('x2')
