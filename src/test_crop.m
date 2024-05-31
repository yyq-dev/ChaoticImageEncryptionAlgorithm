% 裁剪攻击
clc;clear;close all;
P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);
figure(1);imshow(P);
K=[0.1,0.2,0.3,0.4]; % 密钥

C=tpencrypt(P,K);
figure(2);imshow(C);

C1=crop(C,5); % 裁剪比例
figure(3);imshow(C1);

P1=tpdecrypt(C1,K);
figure(4);imshow(P1);