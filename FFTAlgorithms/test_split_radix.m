nu = 5;
if nu<2
    error('nu must be at least 2 for split-radix')
end
N = 2^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
cN = exp(-1j*2*pi/N);
[X_bitrev,n_cadds,n_cmults] = split_radix(x,0,0);
X = decimate_by_2(X_bitrev);

n_radds = 2*(n_cadds+n_cmults);
n_rmults = 4*n_cmults;

diff = fft(x)-X;
energy = sum(abs(diff).^2)
%n_radds/(N*nu) % should be 8/3 ~ 2.66 for large N
%n_rmults/(N*nu) % should be 4/3 ~ 1.33 for large N
