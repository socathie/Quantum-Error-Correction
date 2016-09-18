t = .005;

DDf = [];
DFSf = [];
DFSDDf = [];

for N = 2:6
    QubitBath;
    [~,DDf(N),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,0);
    
    if ~rem(N,2)
        [~,DFSf(N)] = DFS(rho_b,Bx,By,Bz,N,Nb,t);
        [~,DFSDDf(N)] = DFSDD(rho_b,Bx,By,Bz,N,Nb,t,0);
    end
end

DFSf = real(DFSf(2:2:end));
DFSDDf = real(DFSDDf(2:2:end));
DDf = DDf(2:end);

plot(linspace(2,N,N/2),DFSf);
hold on;
plot(linspace(2,N,N-1),DDf);
plot(linspace(2,N,N/2),DFSDDf);