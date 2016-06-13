nu = 5;
N = 3^nu;
n = linspace(0,N-1,N);
x = rand(N,1)+1j*rand(N,1);
n_cadds = 0;
n_cmults = 0;
input = decimate_by_3(x);
output = zeros(N,1);
for w=1:nu
    for d=1:3^(nu-w)
        for m=1:3^w
            if m<=3^(w-1)
                %fprintf('A: w = %i, d = %i, m = %i, a = %i, b = %i, c = %i, (m-1) = %i, 2*(m-1) = %i \n',...
                %        w,d,m,3^w*(d-1)+m,3^w*(d-1)+m+3^(w-1),3^w*(d-1)+m+2*3^(w-1),m-1,2*(m-1))
                output(3^w*(d-1)+m) = input(3^w*(d-1)+m)...
                                      + exp(-1j*2*pi*(m-1)/(3^w))*input(3^w*(d-1)+m+3^(w-1))...
                                      + exp(-1j*2*pi*2*(m-1)/(3^w))*input(3^w*(d-1)+m+2*3^(w-1));
            elseif m<=2*3^(w-1)
                %fprintf('B: w = %i, d = %i, m = %i, a = %i, b = %i, c = %i, (m-1) = %i, 2*(m-1) = %i \n',...
                %        w,d,m,3^w*(d-1)+m-3^(w-1),3^w*(d-1)+m,3^w*(d-1)+m+3^(w-1),m-1,2*(m-1))
                output(3^w*(d-1)+m) = input(3^w*(d-1)+m-3^(w-1))...
                                      + exp(-1j*2*pi*(m-1)/(3^w))*input(3^w*(d-1)+m)...
                                      + exp(-1j*2*pi*2*(m-1)/(3^w))*input(3^w*(d-1)+m+3^(w-1));
            else
                %fprintf('C: w = %i, d = %i, m = %i, a = %i, b = %i, c = %i, (m-1) = %i, 2*(m-1) = %i \n',...
                %        w,d,m,3^w*(d-1)+m-2*3^(w-1),3^w*(d-1)+m-3^(w-1),3^w*(d-1)+m,m-1,2*(m-1))
                output(3^w*(d-1)+m) = input(3^w*(d-1)+m-2*3^(w-1))...
                                      + exp(-1j*2*pi*(m-1)/(3^w))*input(3^w*(d-1)+m-3^(w-1))...
                                      + exp(-1j*2*pi*2*(m-1)/(3^w))*input(3^w*(d-1)+m);               
            end
            n_cmults = n_cmults+2;
            n_cadds = n_cadds+2;
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
n_radds/(N*nu) % should be 8
n_rmults/(N*nu) % should be 8
