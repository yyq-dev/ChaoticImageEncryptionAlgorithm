clc;clear;close all;

P = imread('../pic/v1.png'); % 以v1图像为例
P = rgb2gray(P);

K1=[0.1,0.2,0.3,0.4];
K2=[0.1,0.2,0.3,0.4+10^(-13)]; % 以y1为例

tic;C1=tpencrypt(P,K1);toc
% figure(2);imshow(C1)

tic;C2=tpencrypt(P,K2);toc
% figure(3);imshow(C2)

a = npcruacibaci(C1,C2) % 密钥敏感性