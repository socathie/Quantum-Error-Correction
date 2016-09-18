N = 5;
n = 10; %time intervals
T = linspace(.000005,.00005,n);
QubitBath;
%F = 0*T;
Fs = 0*T;
m = 10; %no. of realizations;

parfor J = 1:n
    t = T(J);
    for I = 1:m
        %F(J) = F(J) + QEC(rho_b,Bx,By,Bz,N,Nb,t);
        Fs(J) = Fs(J) + QECs(rho_b,Bx,By,Bz,N,Nb,t);
    end
end
%F = F/m;
Fs = Fs/m;

%plot(T,F);
%hold on;
plot(T,Fs);