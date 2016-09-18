N = 4;
T = linspace(.001,.01,100);

DDf = [];
CDDf = [];
DFSf = [];

for J = 1:100
    t = T(J);
    QubitBath;
    [~,DDf(J),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,0);
    [~,CDDf(J),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,1);
    [~,DFSf(J)] = DFS(rho_b,Bx,By,Bz,N,Nb,t);
end
DFSf = real(DFSf);

plot(T,DFSf);
hold on;
plot(T,DDf);
plot(T,CDDf);