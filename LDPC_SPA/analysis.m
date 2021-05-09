clear all;
M=csvread('analysis.csv');
M=int16(M');
s=sum(M);
s=(s~=0);
error_rate=sum(s)/length(M);
