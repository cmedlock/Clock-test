function[x] = decimate_by_2(x)

N = length(x);
temp = zeros(N,1);
% decimate once
for m=1:2
    for d=1:N/2
        %fprintf('m=%i,d=%i,n1=%i,n2=%i \n',m,d,d+(m-1)*N/2,m+(d-1)*2)
        temp(d+(m-1)*N/2) = x(m+(d-1)*2);
    end
end
x = temp;
% if each subsequence has N > 2, keep going
if N/2>2
    %fprintf('decimating again \n')
    temp(1:N/2) = decimate_by_2(temp(1:N/2));
    temp(N/2+1:N) = decimate_by_2(temp(N/2+1:N));
end
x = temp;