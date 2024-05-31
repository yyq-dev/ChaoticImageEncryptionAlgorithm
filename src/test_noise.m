% 噪声攻击
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
figure(1);imshow(P);
K=[0.1,0.2,0.3,0.4]; % 密钥

C=tpencrypt(P,K);
figure(2);imshow(C);

C1=imnoise(C,'salt & pepper',0.01); % 椒盐噪声
% C1=imnoise(C,"gaussian",0.0001); % 高斯噪声
figure(3);imshow(C1);

P1=tpdecrypt(C1,K);
figure(4);imshow(P1);