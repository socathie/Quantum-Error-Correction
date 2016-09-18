N = 4;
QubitBath;
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

%free evolution
T = linspace(0,.05,50);
F = [];
for J = 1:50
    t = T(J);
    U = expm(-1i*t*H);
    rho_f = U*rho*U';
    rho_f = TrX23(rho_f,2,[2^N Nb]);
    F(J) = real(sqrt(v'*rho_f*v));
end