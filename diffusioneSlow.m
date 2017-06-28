function F = diffusioneSlow(x,tEnd,u0)
%apply gaussian smoothing to one dimensional signals
%Input variables:
%x = array with spatial discretization (row array with length n+2)
%u0 = initial condition (column row with length n+2)
%tEnd = width smoothing kernel
%Output variables:
%F = value of the smoothed signal

n = length(x);
t = linspace(0, tEnd, 5);
m = length(t)-2;
h = x(2) - x(1);
k = t(2) - t(1);

%set boundary and initial conditions
fIn = ones(1, m+2)*u0(1);
fEnd = ones(1, m+2)*u0(end);
F=u0;


for i=2:m+2
    
    %compute finite difference matrix
    v0 = [0; 2*ones(n-2,1); 0];
    v1 = [0; ones(n-2,1)];
    v2 = [ones(n-2,1); 0];
    %matrix of semimplicit method
    r=(1:n);
    c=(1:n);
    A = sparse(r,c,ones(n,1)+v0*k/((h^2)));
    A = A + sparse(r(1:end-1),c(2:end),-v1*k/((h^2)),n,n);
    A = A + sparse(r(2:end),c(1:end-1),-v2*k/((h^2)),n,n);
    
    %compute next iteration
    F = A\F;
    
    %update values
    F(1) = fIn(i);
    F(end) = fEnd(i);
end