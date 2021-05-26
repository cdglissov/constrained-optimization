%% Damped BFGS
close all
clc
clear 

y0 = [-1,1]';
points = horzcat([0;0],[3;2],[-4.5;3], [-3;-4]);
option1 = {'o-r', 'd-b', '^-g', 's-m'};          
option2 = {'r', 'b', 'g', 'm'};   


PlotMe([-5,5],[-5,5])
for i = 1:length(points)
    [x, stat] = SQP_BFGS_ineq(@ObjFun, @ConFun, points(:,i), y0);
     p(i) = plot(stat.X(1,:),stat.X(2,:),option1{i}, 'linewidth',1.2, ...
        'MarkerFaceColor', option2{i}, 'MarkerSize', 6);
    if i == 1
        table1{i} = table(horzcat((0:(stat.iter))',stat.X'));
    else
        table1{i} = table(stat.X');
    end
end
xlabel('$x_1$','Interpreter','LaTex','FontSize',16);
ylabel('$x_2$','Interpreter','LaTex','FontSize',16);
title("Iteration Sequence, SQP Damped BFGS + Line Search",...
    'Interpreter','LaTex','FontSize',16)
legend([p(1),p(2), p(3), p(4)], ...
    strcat('x_0=(', num2str(points(1,1)),',',num2str(points(2,1)),')'),...
    strcat('x_0=(', num2str(points(1,2)),',',num2str(points(2,2)),')'),...
    strcat('x_0=(', num2str(points(1,3)),',',num2str(points(2,3)),')'),...
    strcat('x_0=(', num2str(points(1,4)),',',num2str(points(2,4)),')'),...
    'Location','SouthEast')
hold off;
%% Line Search + Damped BFGS
close all
clc
clear 

y0 = [-1,1]';
points = horzcat([0;0],[3;2],[-4.5;3], [-3;-4]);
option1 = {'o-r', 'd-b', '^-g', 's-m'};          
option2 = {'r', 'b', 'g', 'm'};   


PlotMe([-5,5],[-5,5])
for i = 1:length(points)
    [x, stat] = SQP_BFGS_LS_ineq(@ObjFun, @ConFun, points(:,i), y0);
     p(i) = plot(stat.X(1,:),stat.X(2,:),option1{i}, 'linewidth',1.2, ...
        'MarkerFaceColor', option2{i}, 'MarkerSize', 6);
    tableStats1{i}=stat;
    if i == 1
        table2{i} = table(horzcat((0:(stat.iter))',stat.X'));
    else
        table2{i} = table(stat.X');
    end
end
xlabel('$x_1$','Interpreter','LaTex','FontSize',16);
ylabel('$x_2$','Interpreter','LaTex','FontSize',16);
title("Iteration Sequence, SQP Damped BFGS + Line Search",...
    'Interpreter','LaTex','FontSize',16)
legend([p(1),p(2), p(3), p(4)], ...
    strcat('x_0=(', num2str(points(1,1)),',',num2str(points(2,1)),')'),...
    strcat('x_0=(', num2str(points(1,2)),',',num2str(points(2,2)),')'),...
    strcat('x_0=(', num2str(points(1,3)),',',num2str(points(2,3)),')'),...
    strcat('x_0=(', num2str(points(1,4)),',',num2str(points(2,4)),')'),...
    'Location','SouthEast')
hold off;

%%
[tableStats1{1}.F(end), tableStats1{1}.iter ,tableStats1{1}.nfun];
[tableStats1{2}.F(end), tableStats1{2}.iter ,tableStats1{2}.nfun];
[tableStats1{3}.F(end), tableStats1{3}.iter ,tableStats1{3}.nfun];
[tableStats1{4}.F(end), tableStats1{4}.iter ,tableStats1{4}.nfun];
%% Trust Region
close all
clc
clear 

y0 = [-1,1]';
points = horzcat([0;0],[3;2],[-4.5;3], [-3;-4]);
option1 = {'o-r', 'd-b', '^-g', 's-m'};          
option2 = {'r', 'b', 'g', 'm'};   

PlotMe([-5,5],[-5,5])
for i = 1:length(points)
    [x, stat] = SQP_TRUST_REG2(@ObjFun, @ConFun, points(:,i), y0);
     p(i) = plot(stat.X(1,:),stat.X(2,:),option1{i}, 'linewidth',1.2, ...
        'MarkerFaceColor', option2{i}, 'MarkerSize', 6);
    tableStats2{i}=stat;
    if i == 1
        table3{i} = table(horzcat((0:(stat.iter))',stat.X'));
    else
        table3{i} = table(stat.X');
    end
end
xlabel('$x_1$','Interpreter','LaTex','FontSize',16);
ylabel('$x_2$','Interpreter','LaTex','FontSize',16);
title("Iteration Sequence, SQP Trust Region",...
    'Interpreter','LaTex','FontSize',16)
legend([p(1),p(2), p(3), p(4)], ...
    strcat('x_0=(', num2str(points(1,1)),',',num2str(points(2,1)),')'),...
    strcat('x_0=(', num2str(points(1,2)),',',num2str(points(2,2)),')'),...
    strcat('x_0=(', num2str(points(1,3)),',',num2str(points(2,3)),')'),...
    strcat('x_0=(', num2str(points(1,4)),',',num2str(points(2,4)),')'),...
    'Location','SouthEast')
hold off;

[tableStats2{1}.F(end), tableStats2{1}.iter ,tableStats2{1}.nfun, tableStats2{1}.ErrL]
[tableStats2{2}.F(end), tableStats2{2}.iter ,tableStats2{2}.nfun, tableStats2{2}.ErrL];
[tableStats2{3}.F(end), tableStats2{3}.iter ,tableStats2{3}.nfun, tableStats2{3}.ErrL];
[tableStats2{4}.F(end), tableStats2{4}.iter ,tableStats2{4}.nfun, tableStats2{4}.ErrL];
