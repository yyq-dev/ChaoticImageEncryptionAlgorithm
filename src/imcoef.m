function r = imcoef(A,N)
    % 图像相关系数计算


    A = double(A); % 将图像矩阵转换为双精度浮点数
    [m, n] = size(A); % 获取图像的行数和列数
    r = zeros(1, 4); % 初始化相关系数向量
    
    % 随机生成像素点的坐标
    x1 = mod(floor(rand(1, N) * 10 ^ 10), m - 1) + 1;
    x2 = mod(floor(rand(1, N) * 10 ^ 10), m) + 1;
    x3 = mod(floor(rand(1, N) * 10 ^ 10), m - 1) + 2;
    y1 = mod(floor(rand(1, N) * 10 ^ 10), n - 1) + 1;
    y2 = mod(floor(rand(1, N) * 10 ^ 10), n) + 1;
    
    % 初始化存储像素值的向量
    u1 = zeros(1, N);
    u2 = zeros(1, N);
    u3 = zeros(1, N);
    u4 = zeros(1, N);
    v1 = zeros(1, N);
    v2 = zeros(1, N);
    v3 = zeros(1, N);
    v4 = zeros(1, N);
    
    % 填充像素值
    for i = 1:N
        u1(i) = A(x1(i), y2(i));
        v1(i) = A(x1(i) + 1, y2(i));
        u2(i) = A(x2(i), y1(i));
        v2(i) = A(x2(i), y1(i) + 1);
        u3(i) = A(x1(i), y1(i));
        v3(i) = A(x1(i) + 1, y1(i) + 1);
        u4(i) = A(x3(i), y1(i));
        v4(i) = A(x3(i) - 1, y1(i) + 1);
    end
    
    % 计算相关系数
    r(1) = mean((u1 - mean(u1)) .* (v1 - mean(v1))) / (std(u1, 1) * std(v1, 1));
    r(2) = mean((u2 - mean(u2)) .* (v2 - mean(v2))) / (std(u2, 1) * std(v2, 1));
    r(3) = mean((u3 - mean(u3)) .* (v3 - mean(v3))) / (std(u3, 1) * std(v3, 1));
    r(4) = mean((u4 - mean(u4)) .* (v4 - mean(v4))) / (std(u4, 1) * std(v4, 1));
    
    % 绘制相关性图形
    figure(101);
    plot(u1, v1, 'b.', 'linewidth', 3, 'markersize', 3);
    axis([0 300 0 300]);
    set(gca, 'XTick', 0:50:300, 'YTick', 0:50:300,'fontsize', 18, 'fontname', 'times new roman');
    set(gca, 'XTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    set(gca, 'YTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    % xlabel('Pixel gray value on location (\itx\rm,\ity\rm)');
    % ylabel('Pixel gray value on location (\itx\rm+1,\ity\rm)');

    figure(102);
    plot(u2, v2, 'b.', 'linewidth', 3, 'markersize', 3);
    axis([0 300 0 300]);
    set(gca, 'XTick', 0:50:300, 'YTick', 0:50:300, 'fontsize', 18, 'fontname', 'times new roman');
    set(gca, 'XTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    set(gca, 'YTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    % xlabel('Pixel gray value on location (\itx\rm,\ity\rm)');
    % ylabel('Pixel gray value on location (\itx\rm,\ity\rm+1)');

    figure(103);
    plot(u3, v3, 'b.', 'linewidth', 3, 'markersize', 3);
    axis([0 300 0 300]);
    set(gca, 'XTick', 0:50:300, 'YTick', 0:50:300, 'fontsize', 18, 'fontname', 'times new roman');
    set(gca, 'XTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    set(gca, 'YTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    % xlabel('Pixel gray value on location (\itx\rm,\ity\rm)');
    % ylabel('Pixel gray value on location (\itx\rm+1,\ity\rm+1)');

    figure(104);
    plot(u4, v4, 'b.', 'linewidth', 3, 'markersize', 3);
    axis([0 300 0 300]);
    set(gca, 'XTick', 0:50:300, 'YTick', 0:50:300, 'fontsize', 18, 'fontname', 'times new roman');
    set(gca, 'XTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    set(gca, 'YTickLabel', {'0', '50', '100', '150', '200', '250', '300'});
    % xlabel('Pixel gray value on location (\itx\rm,\ity\rm)');
    % ylabel('Pixel gray value on location (\itx\rm-1,\ity\rm+1)');

end