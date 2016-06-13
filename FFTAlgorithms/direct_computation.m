nu = 5;
N = 2^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
cN = exp(-1j*2*pi/N);
n_cadds = 0;
n_cmults = 0;
X = zeros(N,1);
for w=0:N-1
    X(w+1) = x(1)*cN^(0*0);
    n_cmults = n_cmults+1;
    for d=1:N-1
        X(w+1) = X(w+1)+x(d+1)*cN^(w*d);
        n_cmults = n_cmults+1;
        n_cadds = n_cadds+1;
    end
end

n_radds = 2*(n_cadds+n_cmults);
n_rmults = 4*n_cmults;

diff = fft(x)-X;
energy = sum(abs(diff).^2)
n_radds/N^2 % should be 4
n_rmults/N^2 % should be 4
