% Driver Exam Question 6
clc; 
clear all;

x_true = [-1.71; 1.59; 1.82; -0.763; -0.763];
x0 = [ -1.8; 1.7; 1.9; -0.8; -0.8];
y0 = [1,1,1]';
objfun = @obj;
confun = @nlpcon;

[x,stat] = NewtonSQP(objfun, confun, x0,y0);

%% 6.2 solve previous with damped BFGS approximation to the Hessian matrix
[x2,stat2] = NewtonSQP_BFGS(objfun, confun, x0,y0);

%% 6.3 solve previous with line search
[x3,stat3] = NewtonSQP_lineSearch(objfun, confun, x0, y0);

%% plots of convergence:
e = linspace(10^-4 , 1);

Err1 = sum(bsxfun(@minus, stat.X, x_true).^2,1);
Err2 = sum(bsxfun(@minus, stat2.X, x_true).^2,1);
Err3 = sum(bsxfun(@minus, stat3.X, x_true).^2,1);

figure(1)
loglog(e , e ,'-','linewidth',1.2);
hold on
loglog(e , e.^2 ,'-','linewidth',1.2);
loglog(Err1(1:(end-1)) , Err1(2:end) , '.-' , 'markersize' , 18,'linewidth',1.2)
loglog(Err2(1:(end-1)) , Err2(2:end) , '.-' , 'markersize' , 18,'linewidth',2)
loglog(Err3(1:(end-1)) , Err3(2:end) , '.-' , 'markersize' , 18,'linewidth',1.2)
hold off
grid on
legend('O(e)' , 'O(e^2)' , 'Local SQP' , 'Local SQP + Damped BFGS',... 
'Local SQP damped BFGS + Line Search','location' , 'southeast' )
xlabel('$\log{e_{k}}$','Interpreter','LaTex','FontSize',18) 
ylabel('$\log{e_{k+1}}$','Interpreter','LaTex','FontSize',18) 
title("Convergence rate of the SQP methods")


%% 
figure(2)
hold on
plot(1:stat.iter+1,stat.Errc,'b.-','linewidth',1.2,'markersize' , 18)
plot(1:stat.iter+1,stat.ErrL,'b*--','linewidth',1.2,'markersize' , 12)
plot(1:stat2.iter+1,stat2.Errc,'r.-','linewidth',1.2,'markersize' , 18)
plot(1:stat2.iter+1,stat2.ErrL,'r*--','linewidth',1.2,'markersize' , 12)
plot(1:stat3.iter+1,stat3.Errc,'g.-','markersize' , 18)
plot(1:stat3.iter+1,stat3.ErrL,'g*-','markersize' , 12)
xlabel('Iterations','FontSize',18)
ylabel('Error','FontSize',18)
legend('SQP\mid\midc(x)\mid\mid_{\infty}','SQP \mid\mid\nabla f(x) - \nabla c(x)^{T}\cdot\lambda\mid\mid_{\infty}',...
    'SQP+BFGS\mid\midc(x)\mid\mid_{\infty}','SQP+BFGS \mid\mid\nabla f(x) - \nabla c(x)^{T}\cdot\lambda\mid\mid_{\infty}',...
    'SQP+BFGS+LS\mid\midc(x)\mid\mid_{\infty}','SQP+BFGS+LS \mid\mid\nabla f(x) - \nabla c(x)^{T}\cdot\lambda\mid\mid_{\infty}')
title("Convergence of the SQP methods" ,'FontSize',18)

%% latex
latex_table1 = latex(vpa(sym(horzcat((1:(stat.iter+1))',stat.X',stat.Y',stat.F')),5));
latex_table2 = latex(vpa(sym(horzcat((1:(stat2.iter+1))',stat2.X',stat2.Y',stat2.F')),5));
latex_table3 = latex(vpa(sym(horzcat((1:(stat3.iter+1))',stat3.X',stat3.Y',stat3.F')),5));

