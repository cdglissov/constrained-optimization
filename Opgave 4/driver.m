%% week6 

%Make Contour plot
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



%% % Markowitz (Test of interior point)

clc;clear;
H = [2.30 0.93 0.62 0.74 -0.23 0;
    0.93 1.40 0.22 0.56 0.26 0;
    0.62 0.22 1.80 0.78 -0.27 0;
    0.74 0.56 0.78 3.40 -0.56 0;
    -0.23 0.26 -0.27 -0.56 2.60 0;
    0 0 0 0 0 0];
mu = [15.1, 12.5, 14.7, 9.02, 17.68, 2]';
e = ones(6,1);
A = [mu,e];
g = 0;
x = ones(size(H,1),1); % 6 variables
y = ones(size(A,2),1);
z = ones(size(x));

C = eye(length(x));
d = zeros(length(x),1);
s = 2.*ones(size(C,1),1);

max1 = 17.68;
min1 = 9.02;
returns = [min1 9.5:0.5:17 max1];
xhats = zeros(5, length(returns));
risk = zeros(length(returns),1);

colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980 ; 0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330]; 

%define max and minimum return
xStar = zeros(6,11);
yStar = zeros(2,11);
zStar = zeros(6,11);
sStar = zeros(6,11);
varStar = zeros(1,11);

colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980 ; 0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330]; 
%define max and minimum return
% Risk free
max1 = 17.68;
min1 = 2;
returns = [min1 3:0.9:16.7 max1];
xhats = zeros(6, length(returns));
risk = zeros(length(returns),1);

for i=1:length(returns)
    b = [returns(i);1];
   [xStar(:,i), yStar(:,i), zStar(:,i), sStar(:,i), info, mu, it] = PDPCIP(H,g,A,C,b,d,x,y,z,s)
    xhats=xStar;
   risk(i) = xStar(:,i)'*H*xStar(:,i);
end

%Risk as function of return
rr = zeros(size(H,1), 1);
figure(1);
plot(returns, risk,'-*b', 'LineWidth', 1)
hold on;
for i= 1:size(xhats,1)
    rr(i) =  plot(A(i,1),H(i,i),'*r','color', colors(i,:),'LineWidth', 2);
end
line([15 15], [0 4], 'Color', 'Black', 'LineWidth', 1.5)
hold off;
xlabel('Return','FontSize',11,'FontSize',14)
ylabel('Risk','FontSize',11,'FontSize',14)
legend(rr, 'Security 1','Security 2','Security 3','Security 4', 'security 5','Security 6','Location','NorthWest')
title('Risk as function of Return', 'FontSize', 14)

% Optimal portfolio as function of return
figure(2);
h2=subplot(4,2,[3 4, 5 ,6,7,8]);
plot(returns, xhats(1,:),'-*b' ,'LineWidth', 1, 'color', colors(1,:))
hold on;
for i=2:6
    plot(returns, xhats(i,:),'-*', 'LineWidth', 1, 'color', colors(i,:))
    set(gca,'FontSize',10)
end
line([15 15], [0 1.07], 'Color', 'Black', 'LineWidth', 1.5)
xlim([min1-0.3, max1+0.3])
ylim([0,1.07])
yticks([0, 0.25, 0.5, 0.75, 1])
xticks([2:16, 17.68])
h1=subplot(4,2,[1,2]);
bar(returns, xhats', 'stacked')
set(gca,'FontSize',10)
legend('Security 1','Security 2','Security 3','Security 4', 'security 5','Security 6','Location','EastOutside')
xlim([min1-0.6, max1+0.6])
ylim([0,1.1])
xticks([2 4.7 7.5 10.2 12.9 15.6, 17.68])
p1=get(h1,'position');
p2=get(h2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1)-0.01 p2(2) p2(3) height],'visible','off');
ylabel('Distribution of Security','visible','on','FontSize',14);
xlabel("Return",'visible','on','FontSize',14)
title("Optimal portfolio as function of return",'visible','on','FontSize',14)




