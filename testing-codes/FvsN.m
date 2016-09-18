t = .005;

DDf = [];
CDDf = [];
DFSf = [];

for N = 2:6
    QubitBath;
    [~,DDf(N),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,0);
    
    [~,CDDf(N),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,1);
    
    if ~rem(N,2)
        [~,DFSf(N)] = DFS(rho_b,Bx,By,Bz,N,Nb,t);
    end
end
DFSf = real(DFSf(2:2:end));
DDf = DDf(2:end);
CDDf = CDDf(2:end);

plot(linspace(2,N,N/2),DFSf);
hold on;
plot(linspace(2,N,N-1),DDf);
plot(linspace(2,N,N-1),CDDf);