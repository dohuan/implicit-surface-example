close all
clear
clc
%%
addpath(genpath('./gpml'))

alpha = (0:0.1:2*pi)';
r = 2;
[on_sur.x,on_sur.y] = pol2cart(alpha,r);

alpha = (0:0.5:2*pi)';
r = 3;
[outer.x,outer.y] = pol2cart(alpha,r);

alpha = (0:0.5:2*pi)';
r = 0.5;
[inner.x,inner.y] = pol2cart(alpha,r);

X = [on_sur.x on_sur.y];
Y = zeros(size(X,1),1);
% X = [X;[inner.x inner.y]];
% Y = [Y;-1*ones(size(inner.x,1),1)];
% X = [X;[outer.x outer.y]];
% Y = [Y;1*ones(size(inner.x,1),1)];

max_x = max(X(:,1));
max_y = max(X(:,2));
min_x = min(X(:,1));
min_y = min(X(:,2));

gridsize = 50;
x_mesh = linspace(max_x,min_x,gridsize);
y_mesh = linspace(max_y,min_y,gridsize);
[S1,S2] = meshgrid(x_mesh,y_mesh); % [middle shortest longest]
S_temp = [S1(:),S2(:)];

covfunc  = @covSEard;
likfunc  = @likGauss;
meanfunc = @meanOne;

hyp.cov(1) = log(0.5); % bandwidth of x
hyp.cov(2) = log(0.5); % bandwidth of y
hyp.cov(3) = log(1);   % sig_f

hyp.lik = log(0.03);

[est, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Y, S_temp);

field = reshape(est,gridsize,gridsize);

figure(1)
surf(S1,S2,field,'EdgeColor','none');
hold on
plot3(on_sur.x,on_sur.y,zeros(size(on_sur.x,1),1),'k','LineWidth',2.5);
set(gca,'FontSize',16);
hold off

figure(2)
pcolor(S1,S2,field);
colorbar
hold on
plot(on_sur.x,on_sur.y,'ko','MarkerSize',10,'LineWidth',2.5);
plot(inner.x,inner.y,'kx','MarkerSize',10,'LineWidth',2.5);
plot(outer.x,outer.y,'kd','MarkerSize',10,'LineWidth',2.5);
hold off
set(gca,'FontSize',16);


