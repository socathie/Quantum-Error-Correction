DFSr = [];
DDr = [];
QECr = [];

for N = 4:100
    DDr(N) = 1;
    QECr(N) = 1-log(1+3*N)/log(2)/N;
    
    if mod(N,2)
        D = nchoosek(N,(N+1)/2);
    else
        D = nchoosek(N,N/2);
    end
    DFSr(N) = log(D)/log(2)/N;
end

plot(linspace(4,N,N-3),DDr(4:end));
hold on
plot(linspace(4,N,N-3),DFSr(4:end));
plot(linspace(4,N,N-3),QECr(4:end));