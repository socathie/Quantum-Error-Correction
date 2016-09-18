%N = 14;
Nb = 2^N;
%t = .005;

x = [0 1;1 0];
y = -[0 1i;-1i 0];
z = [1 0; 0 -1];

rho_b = eye(2^N)/(2^N);
Bx = zeros(2^N);
for j = 1:N
    if j >1
        Id_L = eye(2^(j-1));
    else
        Id_L = 1;
    end
    Id_R = eye(2^(N-j));
    Bx = Bx+kron(Id_L,kron(x,Id_R));
end


By = zeros(2^N);
for j = 1:N
    if j >1
        Id_L = eye(2^(j-1));
    else
        Id_L = 1;
    end
    Id_R = eye(2^(N-j));
    By = By+kron(Id_L,kron(y,Id_R));
end

Bz = zeros(2^N,2^N,N);
for j = 1:N
    if j >1
        Id_L = eye(2^(j-1));
    else
        Id_L = 1;
    end
    Id_R = eye(2^(N-j));
    Bz(:,:,j) = kron(Id_L,kron(z,Id_R));
end