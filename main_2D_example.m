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

hold on
plot(on_sur.x,on_sur.y,'bo','MarkerSize',10);
plot(inner.x,inner.y,'bx','MarkerSize',10);
plot(outer.x,outer.y,'bd','MarkerSize',10);
hold off

X = [on_sur.x on_sur.y];
Y = zeros(size(X,1),1);
X = [X;[inner.x inner.y]];
Y = [Y;-1*ones(size(inner.x,1),1)];
X = [X;[outer.x outer.y]];
Y = [Y;-1*ones(size(inner.x,1),1)];
