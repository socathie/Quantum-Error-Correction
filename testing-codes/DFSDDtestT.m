N = 4;
T = linspace(.001,.01,100);

DDf = [];
DFSf = [];
DFSDDf = [];

for J = 1:100
    t = T(J);
    QubitBath;
    [~,DDf(J),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,0);
    [~,DFSf(J)] = DFS(rho_b,Bx,By,Bz,N,Nb,t);
    [~,DFSDDf(J)] = DFSDD(rho_b,Bx,By,Bz,N,Nb,t,0);
end

DFSf = real(DFSf);
DFSDDf = real(DFSDDf);
%plot(T,DFSf);
hold on;
plot(T,DDf,'x');
plot(T,DFSDDf,'x');

p = polyfit(T,DDf,4); disp(p);
f = polyval(p,T);
plot(T,f);

p = polyfit(T,DFSDDf,4); disp(p);
f = polyval(p,T);
plot(T,f);