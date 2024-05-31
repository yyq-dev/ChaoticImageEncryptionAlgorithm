% 明文与密文直方图及相关系数
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
figure(1);imhist(P) % 明文直方图
K=[0.1,0.2,0.3,0.4]; % 密钥

C=tpencrypt(P,K);
figure(2);imhist(C) % 密文直方图
r1 = imcoef(P,2000) % 明文相关系数 
r2 = imcoef(C,2000) % 密文相关系数