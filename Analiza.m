clc
clear all
%% A weighting
W = csvread('data.csv');
N = W(:,1);
xh = W(:,2);
x = W(:,3);
xs = W(:,4);
e = W(:,5);
eh = W(:,6);
y = W(:,7);