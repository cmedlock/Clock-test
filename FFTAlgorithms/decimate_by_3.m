function[x] = decimate_by_3(x)

N = length(x);
temp = zeros(N,1);
% decimate once
for m=1:3
    for d=1:N/3
        %fprintf('m=%i,d=%i,n1=%i,n2=%i \n',m,d,d+(m-1)*N/3,m+(d-1)*3)
        temp(d+(m-1)*N/3) = x(m+(d-1)*3);
    end
end
x = temp;
% if each subsequence has N > 3, keep going
if N/3>3
    %fprintf('decimating again \n')
    temp(1:N/3) = decimate_by_3(temp(1:N/3));
    temp(N/3+1:2*N/3) = decimate_by_3(temp(N/3+1:2*N/3));
    temp(2*N/3+1:N) = decimate_by_3(temp(2*N/3+1:N));
end
x = temp;