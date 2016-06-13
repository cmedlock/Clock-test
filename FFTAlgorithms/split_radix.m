function[X,n_cadds,n_cmults] = split_radix(x,n_cadds,n_cmults)

N = length(x);
X = zeros(N,1);
[x_split,n_cadds,n_cmults] = split(x,n_cadds,n_cmults);
% even samples of X[k]
if N/2==1
    X(1) = x_split(1);
elseif N/2==2
    X(1) = x_split(1)+x_split(2);
    X(2) = x_split(1)-x_split(2);
    n_cadds = n_cadds+2;
elseif N/2>2
    [X(1:N/2),n_cadds,n_cmults] = split_radix(x_split(1:N/2),n_cadds,n_cmults);
end
% odd samples of X[k]
if N/4==1
    X(3) = x_split(3);
    X(4) = x_split(4);
elseif N/4==2
    X(5) = x_split(5)+x_split(6);
    X(6) = x_split(5)-x_split(6);
    X(7) = x_split(7)+x_split(8);
    X(8) = x_split(7)-x_split(8);
    n_cadds = n_cadds+4;
elseif N/4>2
    [X(N/2+1:3*N/4),n_cadds,n_cmults] = split_radix(x_split(N/2+1:3*N/4),n_cadds,n_cmults);
    [X(3*N/4+1:N),n_cadds,n_cmults] = split_radix(x_split(3*N/4+1:N),n_cadds,n_cmults);
end
    