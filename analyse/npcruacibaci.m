function nu = npcruacibaci(P1, P2)
    % 初始化结果向量
    nu = zeros(1,3);
    % 将输入的图像矩阵转换为double类型
    P1 = double(P1);
    P2 = double(P2);
    % 获取图像大小
    [M, N] = size(P1);
    % 计算NPCR
    D = (P1 ~= P2);
    nu(1) = sum(sum(D)) / (M * N) * 100;
    fprintf('NPCR = %8.4f%%.\n',nu(1));
    % 计算UACI
    nu(2) = sum(sum(abs(P1 - P2))) / (255 * M * N) * 100;
    fprintf('UACI = %8.4f%%.\n',nu(2));
    % 初始化BACI计算变量
    D = abs(P1 - P2);
    m = 0;
    % 遍历图像矩阵计算BACI
    for i = 1:M-1
        for j = 1:N-1
            d = D(i:i+1, j:j+1);
            m = m + (abs(d(1,1) - d(1,2)) + abs(d(1,1) - d(2,1)) + abs(d(1,1) - d(2,2)) + abs(d(1,2) - d(2,1)) + abs(d(1,2) - d(2,2)) + abs(d(2,1) - d(2,2))) / 6 / 255;
        end
    end
    % 计算BACI
    nu(3) = m / ((M - 1) * (N - 1)) * 100;
    fprintf('BACI = %8.4f%%.\n',nu(3));
end
