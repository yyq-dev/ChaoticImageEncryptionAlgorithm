% 信息熵
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
K=[0.1,0.2,0.3,0.4]; % 密钥

C=tpencrypt(P,K);
H1 = entropy(P) % 明文
H2 = entropy(C) % 密文