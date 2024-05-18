function C = tpencrypt(P,K)
%TPENCRYPT 加密图像
%   K = [x0, y0, x1, y1]
% P为输入的明文图像，K为输入密钥，C为输出的密文图像

    [M, N] = size(P); P = double(P); n = 3 * M * N; s = zeros(1, n);
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
    % 序列X
    X = mod(floor((s(1:M * N) + 100) * 10 ^ 10), M * N) + 1; 
    % 序列XT
    XT = mod(floor((s(M * N + 1:2 * M * N) + 100) * 10 ^ 10), M * N) + 1; 
    % 序列SBT
    SBT = mod(floor((s(2*M*N+1:3*M*N)+100)*10^10),256)+1;

    [~, idx] = unique(X); X1 = zeros(1, M * N);
    X1(1 : length(idx)) = X(sort(idx)); 
    X1(length(idx) + 1 : M * N) = setdiff(1:M * N, X1); 
    X = X1;

    [~, idx] = unique(XT); XT1 = zeros(1, M * N);
    XT1(1:length(idx)) = XT(sort(idx)); 
    XT1(length(idx) + 1:M * N) = setdiff(1:M * N, XT1);
    XT = XT1;
    
    % S-Box构建
    keys = [K(3),K(4)];
    S = sboxgen(M,N,keys); % S-Box构建结果

    % 二维图像展开成一维向量后的无重复置乱
    A = P(:);
    for i = 1:floor(M * N / 2)
        t = A(X(i)); A(X(i)) = A(X(M * N - i + 1)); A(X(M * N - i + 1)) = t;
    end

    % 利用生成的S-Box进行扩散处理
    sx = S(:);
    B = zeros(M, N); C = zeros(M, N);
    B0=0;B(1)=bitxor(bitxor(B0,sx(1)),A(1));
    for i=2:M*N
        B(i) = bitxor(bitxor(B(i-1),sx(SBT(i))),A(i));
    end

    C0=0;C(M*N)=bitxor(bitxor(C0,sx(SBT(M*N))),B(M*N));
    for i=M*N-1:-1:1
        C(i) = bitxor(bitxor(C(i+1),sx(SBT(i))),B(i));
    end

    % 无重复置乱
    for i = 1:floor(M * N / 2)
        t = C(XT(i)); C(XT(i)) = C(XT(M * N - i + 1)); C(XT(M * N - i + 1)) = t;
    end

    C=reshape(C,M,N);C=uint8(C);

end

