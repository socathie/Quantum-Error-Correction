function F = QEC(rho_b,Bx,By,Bz,N,Nb,t)
%this encodes 5 physics qubits into 1 encoded qubit using "perfect" code
%rho_f = final system density matrix
%F = fidelity
%t = pulse/gate interval
%rho_b = bath density matrix at time 0
%Bx,By = Nb*Nb matrix
%Bz(:,:,i) = Bz,i => each an Nb*Nb matrix, i = 1 to N
%Nb = dimension of Hilbert space of bath

F = 0;
x = [0 1;1 0];
y = -[0 1i;-1i 0];
z = [1 0; 0 -1];

Sx = zeros(2^N);
Sy = zeros(2^N);
Sz = zeros(2^N,2^N,N);
parfor j = 1:N
    Sz(:,:,j) = zeros(2^N);
end

parfor j = 1:N
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
parfor j = 1:N
    H = H + kron(Sz(:,:,j),Bz(:,:,j));
end

%unitary evolution between gates
U = expm(-1i*t*H);

%total stabilizer generators
S1 = kron(kron(kron(kron(kron(x,z),z),x),eye(2)),eye(Nb));
[V1,D1] = eig(S1);
S2 = kron(kron(kron(kron(kron(eye(2),x),z),z),x),eye(Nb));
[V2,D2] = eig(S2);
S3 = kron(kron(kron(kron(kron(x,eye(2)),x),z),z),eye(Nb));
[V3,D3] = eig(S3);
S4 = kron(kron(kron(kron(kron(z,x),eye(2)),x),z),eye(Nb));
[V4,D4] = eig(S4);

%single qubit basis state
u = [1;0];
d = [0;1];

%encode into logical zero
v = kron(kron(kron(kron(u,u),u),u),u);
v = v + kron(kron(kron(kron(d,u),u),d),u);
v = v + kron(kron(kron(kron(u,d),u),u),d);
v = v + kron(kron(kron(kron(d,u),d),u),u);
v = v + kron(kron(kron(kron(u,d),u),d),u);
v = v - kron(kron(kron(kron(d,d),u),d),d);
v = v - kron(kron(kron(kron(u,u),d),d),u);
v = v - kron(kron(kron(kron(d,d),u),u),u);
v = v - kron(kron(kron(kron(d,d),d),u),d);
v = v - kron(kron(kron(kron(u,u),u),d),d);
v = v - kron(kron(kron(kron(d,d),d),d),u);
v = v - kron(kron(kron(kron(u,d),d),d),d);
v = v - kron(kron(kron(kron(d,u),u),u),d);
v = v - kron(kron(kron(kron(u,d),d),u),u);
v = v - kron(kron(kron(kron(d,u),d),d),d);
v = v + kron(kron(kron(kron(u,u),d),u),d);
v = v/4;

rho_s = v*v';
rho = kron(rho_s,rho_b);

%evolution
rho_f = U*rho*U';

%1st stablilzer measurement
parfor i = 1:(2^N*Nb)
    proj = V1(:,i)*V1(:,i)';
    p(i) = trace(rho_f*proj);
end
p = real(p);
pr = rand();
temp = 0;
j = 1;
while pr > temp
    temp = temp + p(j);
    j = j+1;
end
s1 = D1(j-1,j-1);
proj = V1(:,j-1)*V1(:,j-1)';
rho_f = proj*rho_f*proj/p(j-1);

%evolution
rho_f = U*rho_f*U';

%2nd stablilzer measurement
parfor i = 1:(2^N*Nb)
    proj = V2(:,i)*V2(:,i)';
    p(i) = trace(rho_f*proj);
end
p = real(p);
pr = rand();
temp = 0;
j = 1;
while pr > temp
    temp = temp + p(j);
    j = j+1;
end
s2 = D2(j-1,j-1);
proj = V2(:,j-1)*V2(:,j-1)';
rho_f = proj*rho_f*proj/p(j-1);

%evolution
rho_f = U*rho_f*U';

%3rd stablilzer measurement
parfor i = 1:(2^N*Nb)
    proj = V3(:,i)*V3(:,i)';
    p(i) = trace(rho_f*proj);
end
p = real(p);
pr = rand();
temp = 0;
j = 1;
while pr > temp
    temp = temp + p(j);
    j = j+1;
end
s3 = D3(j-1,j-1);
proj = V3(:,j-1)*V3(:,j-1)';
rho_f = proj*rho_f*proj/p(j-1);

%evolution
rho_f = U*rho_f*U';

%4th stablilzer measurement
parfor i = 1:(2^N*Nb)
    proj = V4(:,i)*V4(:,i)';
    p(i) = trace(rho_f*proj);
end
p = real(p);
pr = rand();
temp = 0;
j = 1;
while pr > temp
    temp = temp + p(j);
    j = j+1;
end
s4 = D4(j-1,j-1);
proj = V4(:,j-1)*V4(:,j-1)';
rho_f = proj*rho_f*proj/p(j-1);

%syndrome identification
s = [sign(s1), sign(s2), sign(s3), sign(s4)];

if isequal(s,[1 1 1 1]) %"no" error
    corr = eye(2^N*Nb);
elseif isequal(s,[1 1 1 -1]) %X error on 1st qubit
    corr = kron(x,eye(2^N*Nb/2));
elseif isequal(s,[-1 1 -1 -1]) %Y error on 1st qubit
    corr = kron(y,eye(2^N*Nb/2));
elseif isequal(s,[-1 1 -1 1]) %Z error on 1st qubit
    corr = kron(z,eye(2^N*Nb/2));
elseif isequal(s,[-1 1 1 1]) %X error on 2nd qubit
    corr = kron(kron(eye(2),x),eye(2^N*Nb/4));
elseif isequal(s,[-1 -1 1 -1]) %Y error on 2nd qubit
    corr = kron(kron(eye(2),y),eye(2^N*Nb/4));
elseif isequal(s,[1 -1 1 -1]) %Z error on 2nd qubit
    corr = kron(kron(eye(2),z),eye(2^N*Nb/4));
elseif isequal(s,[-1 -1 1 1]) %X error on 3rd qubit
    corr = kron(kron(eye(4),x),eye(2^N*Nb/8));
elseif isequal(s,[-1 -1 -1 1]) %Y error on 3rd qubit
    corr = kron(kron(eye(4),y),eye(2^N*Nb/8));
elseif isequal(s,[1 1 -1 1]) %Z error on 3rd qubit
    corr = kron(kron(eye(4),z),eye(2^N*Nb/16));
elseif isequal(s,[1 -1 -1 1]) %X error on 4th qubit
    corr = kron(kron(eye(8),x),eye(2^N*Nb/16));
elseif isequal(s,[-1 -1 -1 -1]) %Y error on 4th qubit
    corr = kron(kron(eye(8),y),eye(2^N*Nb/16));
elseif isequal(s,[-1 1 1 -1]) %Z error on 4th qubit
    corr = kron(kron(eye(8),z),eye(2^N*Nb/16));
elseif isequal(s,[1 1 -1 -1]) %X or Z error on 5th qubit
    corr = kron(kron(eye(16),x),eye(2^N*Nb/32));
elseif isequal(s,[1 -1 -1 -1]) %Y error on 5th qubit
    corr = kron(kron(eye(16),y),eye(2^N*Nb/32));
else
    disp('fail');
    return
end

%find final system density matrix
rho_f = U*rho_f*U';
%correction
rho_f = corr*rho_f*corr;
%take partial trace
rho_f = TrX23(rho_f,2,[2^N Nb]);
%find fidelity
F = real(sqrt(v'*rho_f*v));