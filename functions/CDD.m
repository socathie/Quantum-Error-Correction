function [rho_f,F,T] = CDD(rho_b,Bx,By,Bz,t,N,Nb,m)
%rho_f = system density matrix at time T
%F = fidelity
%T = total time (in units of t)
%rho_b = bath density matrix at time 0
%Bx,By = Nb*Nb matrix
%Bz(:,:,i) = Bz,i => each an Nb*Nb matrix, i = 1 to N
%t = tau (pulse interval)
%N = system size
%Nb = dimension of Hilbert space of bath
%m = concatenation level (m = 0 => ordinary DD)

%define single qubit Pauli matrices
x = [0 1;1 0];
y = -[0 1i;-1i 0];
z = [1 0; 0 -1];

%initial state = tensor product of all |0> = (1,0) states
v = 1;
for j = 1:N
    v = kron([1;0],v);
end
rho_i = v*v';

%total initial density matrix
rho = kron(rho_i,rho_b);

Sx = zeros(2^N);
Sy = zeros(2^N);
Sz = zeros(2^N,2^N,N);
for j = 1:N
    Sz(:,:,j) = zeros(2^N);
end

for j = 1:N
    if j>1
        Id_L = eye(2^(j-1));
    else
        Id_L = 1;
    end
    Id_R = eye(2^(N-j));
    Sx = Sx + kron(Id_L,kron(x,Id_R));
    Sy = Sy + kron(Id_L,kron(y,Id_R));
    Sz(:,:,j) = kron(Id_L,kron(z,Id_R));
end

%compute the SB Hamiltonian
H = kron(Sx,Bx) + kron(Sy,By);
for j = 1:N
    H = H + kron(Sz(:,:,j),Bz(:,:,j));
end

%free evolution for each pulse interval
U = expm(-1i*t*H);

%collective Z pulses
Zp = eye(2^N);
for j = 1:N
    Zp = Zp*Sz(:,:,j);
end
Zp = kron(Zp,eye(Nb));

%collective X pulses
Xp = eye(2^N);
for j = 1:N
    if j >1
        Id_L = eye(2^(j-1));
    else
        Id_L = 1;
    end
    Id_R = eye(2^(N-j));
    Xp = Xp*kron(Id_L,kron(x,Id_R));
end
Xp = kron(Xp,eye(Nb));
T = t;

%apply pulses in each level
for k = 1:m+1
    T = 4*T;
    U = Zp*U*Zp*U;
    U = Xp*U*Xp*U;
end

%find final system density matrix
rho_f = U*rho*U';
%take partial trace
rho_f = TrX23(rho_f,2,[2^N Nb]);

%find fidelity
F = sqrt(v'*rho_f*v);
