clc;clear;close all;
% 绘制y-c1的分岔图，其他坐标类似

% 参数设置
a1 = 0.9;
b1 = 0.7;
a2 = 1.2;
b2 = 1.5;
c1 = 0.6;
c2 = 0.8;

% 分岔图参数设置
y0 = 0.2;
x0 = 0.3;
iterations = 1000; % 总迭代次数
transients = 500; % 跳过前500次迭代，以去除瞬态效应
c1_range = linspace(0, 10, 1000); % c1参数的变化范围

% 初始化数组
x = zeros(iterations, 1);
y = zeros(iterations, 1);
x(1) = x0;
y(1) = y0;

% 准备绘图
figure;
hold on;

% 分岔图计算
for c1 = c1_range
    % 迭代计算
    for n = 1:iterations-1
        x(n + 1) = mod(exp(a1) * y(n) + exp(b1) * (y(n) + c1)^2, 1);
        y(n + 1) = mod(exp(a2) * x(n) + exp(b2) * (x(n) + c2)^3, 1);
    end
    
    % 去除瞬态效应
    x_trans = y(transients:end);
    
    % 绘制分岔图
    plot(c1 * ones(length(x_trans), 1), x_trans, '.k', 'MarkerSize', 0.5);
end

hold off;
set(gca,'fontsize',30,'FontName','times new roman');
xlabel('c_1','FontName','Times New Roman','FontSize',30);
ylabel('y','FontName','Times New Roman','FontSize',30);

