function P = tpdecrypt(C,K)
%TPDECRYPT 解密图像
%   
% C为输入的密文图像，K为输入密钥，P为输出的明文图像
    
    [M, N] = size(C); C = double(C); n = 3 * M * N; s = zeros(1, n);
    h = 0.002; t = 1000;
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
    X = mod(floor((s(1:M * N) + 100) * 10^10), M * N) + 1;
    XT = mod(floor((s(M * N + 1:2 * M * N) + 100) * 10 ^ 10), M * N) + 1; 
    [~, idx] = unique(X); X1 = zeros(1, M * N);
    X1(1:length(idx)) = X(sort(idx)); 
    X1(length(idx) + 1:M * N) = setdiff(1:M * N, X1);
    X = X1;

    [~, idx] = unique(XT); XT1 = zeros(1, M * N);
    XT1(1:length(idx)) = XT(sort(idx)); 
    XT1(length(idx) + 1:M * N) = setdiff(1:M * N, XT1);
    XT = XT1;

    % S-Box构建
    list(1:256) = s(1:256);list1(1:128) = s(257:257+127);
    list2(1:128) = s(257+127:257+127*2);
    [~,L]=sort(list);[~,L1]=sort(list1);[~,L2]=sort(list2);
    L = mod(L,256);L1 = mod(L1,128);L2 = mod(L2,128);
    L1 = L1 + 128;L2 = L2 + 128;
    S1 = reshape(L,16,16); a1 = S1';

    s1 = a1(:);
    for i=1:128
        t = s1(i);s1(i) = s1(L1(i));s1(L1(i)) = t;
    end

    S2 = reshape(s1,16,16); a2 = S2';
    s2 = a2(:);
    for i=1:128
        t = s2(i);s2(i) = s2(L2(i));s2(L2(i)) = t;
    end

    S = reshape(s2,16,16); % S-Box构建结果

    C = C(:);  % 置乱的逆
    for i = 1:floor(M * N / 2)
        t = C(XT(i)); C(XT(i)) = C(XT(M * N - i + 1)); C(XT(M * N - i + 1)) = t;
    end

    sx = S(:);  % 扩散的逆
    A = C;  
    D = zeros(M, N); E = zeros(M, N);
    DO = 0; 
    D(M*N) = bitxor(bitxor(DO, sx(mod(M*N,256)+1)), A(M*N));
    for i = M*N-1:-1:1
        D(i) = bitxor(bitxor(C(i+1), sx(mod(i,256)+1)), A(i));
    end
    E0 = 0; 
    E(1) = bitxor(bitxor(E0, sx(1)), D(1));
    for i = 2:M*N
        E(i) = bitxor(bitxor(D(i-1), sx(mod(i,256)+1)), D(i));
    end
    
    % 置乱的逆
    for i = 1:floor(M * N / 2)
        t = E(X(i)); E(X(i)) = E(X(M * N - i + 1)); E(X(M * N - i + 1)) = t;
    end
    P = reshape(E, M, N); P = uint8(P);

end

