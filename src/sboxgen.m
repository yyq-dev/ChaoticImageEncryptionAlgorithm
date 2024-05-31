function S = sboxgen(M,N,K)
% 基于混沌序列的S-box生成
    
    n = M * N; s = zeros(1, n); % 混沌序列生成
    h = 0.002; t = 3000;
    a1=0.9;b1=0.7;c1=0.6;
    a2=1.2;b2=1.5;c2=0.8;
    x0 = K(1); y0 = K(2); 
    for i = 1:n + t
        x1=mod(exp(a1)*y0+exp(b1)*(y0+c1)^2,1);
        y1=mod(exp(a2)*x0+exp(b2)*(x0+c2)^3,1);
        x0 = x1; y0 = y1;
        if i > t
            s(i - t) = x1;
            if mod((i - t), 3000) == 0
                x0 = x0 + h * sin(y0);
            end
        end
    end

    % S-Box构建
    list(1:256) = s(1:256);list1(1:128) = s(257:257+127);
    list2(1:128) = s(257+127:257+127*2);
    [~,L]=sort(list);[~,L1]=sort(list1);[~,L2]=sort(list2);
    L = mod(L,256);L1 = mod(L1,128);L2 = mod(L2,128);
    L1 = L1 + 128;L2 = L2 + 128;
    S1 = reshape(L,16,16); a1 = S1';

    s1 = a1(:); % 进行置乱
    for i=1:128
        t = s1(i);s1(i) = s1(L1(i));s1(L1(i)) = t;
    end

    S2 = reshape(s1,16,16); a2 = S2';
    s2 = a2(:);
    for i=1:128
        t = s2(i);s2(i) = s2(L2(i));s2(L2(i)) = t;
    end
    
    S = reshape(s2,16,16); % S-Box构建结果
end