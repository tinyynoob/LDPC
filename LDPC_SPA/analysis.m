clear all;
M=int16(readmatrix('analysis.csv'));
M=M';
s=logical(sum(M));
error_rate=sum(s)/length(s);
