nu = 2;
N = 3^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
n_cadds = 0;
n_cmults = 0;
x_bitrev = decimate_by_3(x);
