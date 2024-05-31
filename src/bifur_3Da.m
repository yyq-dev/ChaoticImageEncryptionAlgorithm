clc;clear;close all;

% 参数设置
b1 = 0.7;
c1 = 0.6;
b2 = 1.5;
c2 = 0.8;

% 分岔图参数设置
y0 = 0.2;
x0 = 0.3;
iterations = 1000; % 总迭代次数
transients = 500; % 跳过前500次迭代，以去除瞬态效应
a1_range = linspace(0, 10, 100); % a1参数的变化范围
a2_range = linspace(0, 10, 100); % a2参数的变化范围

% 初始化数组
x = zeros(iterations, 1);
y = zeros(iterations, 1);
x(1) = x0;
y(1) = y0;

% 预分配结果矩阵
X = zeros(length(a1_range), length(a2_range), iterations - transients);

% 分岔图计算
for i = 1:length(a1_range)
    a1 = a1_range(i);
    for j = 1:length(a2_range)
        a2 = a2_range(j);
        
        % 迭代计算
        for n = 1:iterations-1
            x(n + 1) = mod(exp(a1) * y(n) + exp(b1) * (y(n) + c1)^2, 1);
            y(n + 1) = mod(exp(a2) * x(n) + exp(b2) * (x(n) + c2)^3, 1);
        end
        
        % 去除瞬态效应并存储结果
        X(i, j, :) = x(transients+1:end);
    end
end

% 生成三维图
figure;
hold on;

% 绘制三维分岔图
for k = 1:iterations - transients
    surf(a1_range, a2_range, squeeze(X(:,:,k)), 'EdgeColor', 'none', 'FaceAlpha', 0.1);
end

hold off;
xlabel('a_1');
ylabel('a_2');
zlabel('x_n');

view(3);
