function F = PMconTV(x,t,s0,met,soglia,sigma,delta)
%function to apply the modified PM filter to one dimensional signals
%Input variables:
%x = array with spatial discretization (row array with length n+2)
%t = array with temporal discr tization (row array with length m+2)
%s0 = initial condition (column row with length n+2)
%met = choice of the edge stopping function
%soglia = edge-stopping threshold
%sigma = width smoothing kernel
%delta = width of the interval for the computation of TV e LV
%Output variables:
%F = value of the smoothed function

n = length(x) - 2;
m = length(t) - 2;
h = x(2) - x(1);
k = t(2) - t(1);
eps = 0.001;

%set boundary conditions
fIn = ones(1, m+2)*s0(1);
fEnd = ones(1, m+2)*s0(end);

%set intial and boundary conditions
F = s0;
if (sigma > 0)
    Fnew = diffusioneSlow(x, sigma, F);
    F = Fnew;
end
Festesa = [F; F(end)];

%compute gradient. TV and LV
deltSig = diff(Festesa);
absdeltSig = abs(deltSig);
filtered1 = filter(ones(1, delta), 1, flip(absdeltSig));
filtered2 = abs(filter(ones(1, delta), 1, flip(deltSig)));
filt1 = flip(filtered1);
filt2 = flip(filtered2);

%evaluate edge-stopping function
pGrad = metodo1(0 + (filt2)./abs(eps + filt1), met, soglia);

for i = 2:m+2
    
    %compute finite difference matrix
    f1 = filter([-1 -2 -1], 1, pGrad);
    f2 = filter([1 1], 1, pGrad);
    v0 = [0; f1(3:end); 0];
    v1 = [0; f2(3:end)];
    v2 = [f2(2:end-1); 0];
    %matrix of semimplicit method
    r = (1:n+2);
    c = (1:n+2);
    A = sparse(r, c, ones(n+2, 1)-v0*k / (2*(h^2)));
    A = A + sparse(r(1:end-1), c(2:end), -v1*k / (2*(h^2)), n+2, n+2);
    A = A + sparse(r(2:end), c(1:end-1), -v2*k / (2*(h^2)), n+2, n+2);

    %compute next iteration
    F = A\F;
    
    %update values
    if ( sigma > 0)
        Fnew = diffusioneSlow(x, sigma, F);
        F = Fnew;
    end
    Festesa = [F; F(end)];
    deltSig = diff(Festesa);
    absdeltSig = abs(deltSig);
    filtered1 = filter(ones(1, delta), 1, flip(absdeltSig));
    filtered2 = abs(filter(ones(1, delta), 1, flip(deltSig)));
    filt1 = flip(filtered1);
    filt2 = flip(filtered2);
    %evaluate edge-stopping function
    pGrad = metodo1((filt2)./abs(eps + filt1), met, soglia);
    
    F(1) = fIn(i);
    F(end) = fEnd(i);  
end