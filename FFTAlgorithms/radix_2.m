nu = 5;
N = 2^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
n_cadds = 0;
n_cmults = 0;
input = decimate_by_2(x);
output = zeros(N,1);
for w=1:nu
    for d=1:2^(nu-w)
        for m=1:2^w
            if m<=2^(w-1)
                %fprintf('A: w = %i, d = %i, m = %i, a = %i, b = %i \n',w,d,m,2^w*(d-1)+m,2^w*(d-1)+m+2^(w-1))
                output(2^w*(d-1)+m) = input(2^w*(d-1)+m)+exp(-1j*2*pi*(m-1)/(2^w))*input(2^w*(d-1)+m+2^(w-1));
            else
                %fprintf('B: w = %i, d = %i, m = %i, a = %i, b = %i \n',w,d,m,2^w*(d-1)+m,2^w*(d-1)+m-2^(w-1))
                output(2^w*(d-1)+m) = exp(-1j*2*pi*(m-1)/(2^w))*input(2^w*(d-1)+m)+input(2^w*(d-1)+m-2^(w-1));
            end
            n_cmults = n_cmults+1;
            n_cadds = n_cadds+1;
        end
    end
    input = output;
    fprintf('\n')
end
X = output;

n_radds = 2*(n_cadds+n_cmults);
n_rmults = 4*n_cmults;

diff = fft(x)-X;
energy = sum(abs(diff).^2)
%n_radds/(N*nu) % should be 4
%n_rmults/(N*nu) % should be 4
