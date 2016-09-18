function [rho_f,F] = DFS(rho_b,Bx,By,Bz,N,Nb,T)
%this encodes N physical qubits into N!/[(N/2)!(N/2+1)!] logical qubits
%requires N even, and only works for initial state being logical zero
%rho_f = system density matrix at time T
%F = fidelity
%T = runtime
%rho_b = bath density matrix at time 0
%Bx,By = Nb*Nb matrix
%Bz(:,:,i) = Bz,i => each an Nb*Nb matrix, i = 1 to N
%N = system size
%Nb = dimension of Hilbert space of bath

%single qubit basis state
u = [1;0];
d = [0;1];

%define singlet
s = (kron(u,d)-kron(d,u))/sqrt(2);
v = 1;

%initial system state
for j = 1:N/2
    v = kron(v,s);
end

rho_s = v*v';

%define single qubit Pauli matrices
x = [0 1;1 0];
y = -[0 1i;-1i 0];
z = [1 0; 0 -1];

%total initial density matrix
rho = kron(rho_s,rho_b);

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

%unitary evolution
U = expm(-1i*T*H);

%find final system density matrix
rho_f = U*rho*U';
%take partial trace
rho_f = TrX23(rho_f,2,[2^N Nb]);

%find fidelity
F = sqrt(v'*rho_f*v);
