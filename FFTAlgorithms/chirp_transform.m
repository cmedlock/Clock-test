nu = 5;
N = 2^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
cN = exp(-1j*2*pi/N);
n_cadds = 0;
n_cmults = 0;
X = zeros(N,1);

w_0 = 0;
delta_w = 2*pi/N;
for k=0:N-1
    conv_sum = 0;
    for m=0:N-1
        x_0 = x(m+1)*exp(-1j*w_0*m)*exp(-1j*delta_w/2*m^2);
        n_cmults = n_cmults+1;
        conv_sum = conv_sum+x_0*exp(1j*delta_w/2*(k-m)^2);
        if m>0
            n_cadds = n_cadds+1;
        end
    end
    X(k+1) = conv_sum*exp(-1j*delta_w/2*k^2);
end

n_radds = 2*(n_cadds+n_cmults);
n_rmults = 4*n_cmults;

diff = fft(x)-X;
energy = sum(abs(diff).^2)
n_radds/N^2 % should be 4 for large N
n_rmults/N^2 % should be 4 for large N
