% 加解密时间
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
figure(1);imshow(P)
K=[0.1,0.2,0.3,0.4]; % 密钥

tic;C=tpencrypt(P,K);toc
figure(2);imshow(C)

tic;P1=tpdecrypt(C,K);toc
figure(3);imshow(P1)