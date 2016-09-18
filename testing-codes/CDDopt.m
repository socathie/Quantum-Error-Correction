N = 4;
t = .005;
QubitBath;
F = [];

for j = 1:10
    [~,F(j),~] = CDD(rho_b,Bx,By,Bz,t,N,Nb,j);
end

%plot(F);

[~,I] = max(F);
disp(I);