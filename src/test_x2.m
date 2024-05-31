% 图像直方图x^2检验
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
K=[0.1,0.2,0.3,0.4]; % 密钥

C=tpencrypt(P,K);

[M,N] = size(P); g = M*N/256;
P = double(P);C = double(C);
fp1 = hist(P(:),256);chai1 = sum((fp1 - g).^2)/g  % 明文
fc1 = hist(C(:),256);chaic1 = sum((fc1 - g).^2)/g % 密文