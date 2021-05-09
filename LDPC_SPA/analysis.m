clear all;
M = int8(readmatrix('analysis.csv'));
u = size(M);
totalBitNum = u(1)*(u(2)-1);    %subtract the extra column

M = M';
s = sum(M);
bit_error_rate = sum(s) / totalBitNum;
frame_error_rate = sum(logical(s)) / length(s);
