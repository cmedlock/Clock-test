nu = 10;
N = 2^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
cN = exp(-1j*2*pi/N);
n_radds = 0;
n_rmults = 0;
X = zeros(N,1);
for k=0:N/2
    w_nminusone = 0;
    w_nminustwo = 0;
    w_n = 0;
    for n=0:N
        if n==N
            w_n = 0+2*cos(2*pi*k/N)*w_nminusone-w_nminustwo;
            n_rmults = n_rmults+2;
            n_radds = n_radds+4;
            break
        else
            w_n = x(n+1)+2*cos(2*pi*k/N)*w_nminusone-w_nminustwo;
            n_rmults = n_rmults+2;
            n_radds = n_radds+4;
        end
        w_nminustwo = w_nminusone;
        w_nminusone = w_n;
    end
    X(k+1) = w_n-cN^k*w_nminusone;
    n_rmults = n_rmults+4;
    n_radds = n_radds+4;
    if k>0 & k<N/2
        X(N-k+1) = w_n-cN^(N-k)*w_nminusone;
        n_rmults = n_rmults+4;
        n_radds = n_radds+4;
    end
end

diff = fft(x)-X;
energy = sum(abs(diff).^2)
n_radds/N^2 % should be 2 for large N
n_rmults/N^2 % should be 1 for large N
