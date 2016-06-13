function[x,n_cadds,n_cmults] = split(x,n_cadds,n_cmults)

N = length(x);
cN = exp(-1j*2*pi/N);
temp = zeros(N,1);
for m=1:N
    if m<=N/2
        temp(m) = x(m)+x(m+N/2);
        n_cadds = n_cadds+1;
    else
        temp(m) = x(m-N/2)-x(m);
        n_cadds = n_cadds+1;
    end
end
temp2 = temp;
for m=N/2+1:N
    if m<=3*N/4
        temp(m) = (temp2(m)-1j*temp2(m+N/4))*cN^(m-N/2-1);
        n_cadds = n_cadds+1;
        n_cmults = n_cmults+1;
    else
        temp(m) = (temp2(m-N/4)+1j*temp2(m))*cN^(3*(m-3*N/4-1));
        n_cadds = n_cadds+1;
        n_cmults = n_cmults+1;        
    end
end
x = temp;
