%% Opgave 3.3
%Specify return
R=10;
[H,g,A,b,C,d] = ConstructMarkowitz(R,0);
%Solve QP using quadprog
xhat = quadprog(H,[],-C,d,A',b);

%find variance
xhat'*H*xhat;

%% Opgave 3.4
colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980 ; 0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330]; 
%define max and minimum return
max1 = 17.68;
min1 = 9.02;
returns = [min1 9.5:0.5:17 max1];
xhats = zeros(5, length(returns));
risk = zeros(length(returns),1);

for i=1:length(returns)
    [H,g,A,b,C,d] = ConstructMarkowitz(returns(i),0);
    xhat = quadprog(H,[],-C,d,A',b);
    xhats(:,i)=xhat;
    risk(i) = xhat'*H*xhat;
end

%Risk as function of return
figure(1);
plot(returns, risk,'-*b', 'LineWidth', 1)
xlabel('Return','FontSize',11)
ylabel('Risk','FontSize',11)
title('Risk as function of Return', 'FontSize', 10)

% Optimal portfolio as function of return
figure(2);
h1=subplot(4,2,[1,2]);
bar(returns, xhats', 'stacked')
set(gca,'FontSize',10)
legend('Security 1','Security 2','Security 3','Security 4', 'security 5','Location','EastOutside')
xlim([min1-0.3, max1+0.3])
ylim([0,1.1])
xticks([9.02, 10:16, 17.68])
h2=subplot(4,2,[3 4, 5 ,6,7,8]);
plot(returns, xhats(1,:),'-*b' ,'LineWidth', 1)
hold on;
for i=2:5
    plot(returns, xhats(i,:),'-*', 'LineWidth', 1,'color', colors(i,:))
    set(gca,'FontSize',10)
end
hold off;
xlim([min1-0.3, max1+0.2])
yticks([0, 0.25, 0.5, 0.75, 1])
xticks([9.02, 10:17, 17.68])
p1=get(h1,'position');
p2=get(h2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1)-0.01 p2(2) p2(3) height],'visible','off');
ylabel('Distribution of Security','visible','on');
xlabel("Return",'visible','on')
title("Optimal portfolio as function of return",'visible','on')

%%  Risk free
max1 = 17.68;
min1 = 2;
returns = [min1 3:0.9:16.7 max1];
xhats = zeros(6, length(returns));
risk = zeros(length(returns),1);

for i=1:length(returns)
    [H,g,A,b,C,d] = ConstructMarkowitz(returns(i),1);
    xhat = quadprog(H,[],-C,d,A',b);
    xhats(:,i)=xhat;
    risk(i) = xhat'*H*xhat;
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
xlabel('Return','FontSize',11)
ylabel('Risk','FontSize',11)
legend(rr, 'Security 1','Security 2','Security 3','Security 4', 'security 5','Security 6','Location','NorthWest')
title('Risk as function of Return', 'FontSize', 10)

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
ylabel('Distribution of Security','visible','on');
xlabel("Return",'visible','on')
title("Optimal portfolio as function of return",'visible','on')

%%
r=15;
[H,g,A,b,C,d] = ConstructMarkowitz(r,1);
xhat = quadprog(H,[],-C,d,A',b);
risk = xhat'*H*xhat;

