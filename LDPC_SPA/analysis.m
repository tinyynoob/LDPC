clear all;
M=csvread('analysis.csv');
s=zeros(length(M),1);
for i=1:length(M)
    s(i)=sum(M(i,:));
end
s=(s~=0);
error_rate=sum(s)/length(M);
