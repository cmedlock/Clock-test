function[x,n_cadds,n_cmults] = split_v2(x,n_cadds,n_cmults)

N = length(x);
cN = exp(-1j*2*pi/N);
temp = zeros(N,1);
for m=1:N
    if m<=N/2
        temp(m) = x(m)+x(m+N/2);
    elseif m<=3*N/4
        temp(m) = ((x(m-N/2)-x(m))-1j*(x(m-N/4)-x(m+N/4)))*cN^(m-N/2-1);
    else
        fprintf('m=%i,3*(m-3*N/4-1)=%i \n',m,3*(m-3*N/4-1))
        temp(m) = ((x(m-3*N/4)-x(m-N/4))+1j*(x(m-N/2)-x(m)))*cN^(3*(m-3*N/4-1));
    end
end
x = temp;
