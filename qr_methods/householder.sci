function [U, R] = qr_house_v1(A)
    [m, n] = size(A);
    k = min(m,n); //define sobre quem iterar
    U = zeros(m, n);
    R = A;

    for j = 1:k //se j > m ou j > n para
        x = R(j:m, j); //elementos abaixo da diagonal

        e = zeros(length(x), 1);
        e(1) = norm(x);

        //para evitar erros númericos
        if x(1) > 0
            u = x + e;
        else
            u = x - e;
        end

        u = u / norm(u); //normaliza u

        //transformação de Householder
        R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n));
        U(j:m, j) = u;
    end
endfunction

function [U, R] = qr_house_v2(A)
    [m, n] = size(A);
    k = min(m-1, n); //define sobre quem iterar (m-1)
    U = zeros(m, k);
    R = A;

    for j = 1:k //se j > m-1 ou j > n para
        x = R(j:m, j); //elementos abaixo da diagonal

        e = zeros(length(x), 1);
        e(1) = norm(x);

        //para evitar erros númericos
        if x(1) > 0
            u = x + e;
        else
            u = x - e;
        end
        
        u = u / norm(u); //normaliza u

        //transformação de Householder
        R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n));
        U(j:m, j) = u;
    end
endfunction

function Q = constroi_Q_House(U)
    [m, k] = size(U);
    Q = eye(m, m);

    //percorre do fim para o começo da matriz
    for j = k:-1:1
        u = U(:, j); //seleciona a coluna

        //atualiza Q aplicando a transformação de Householder
        Q = Q - 2 * u * (u' * Q);
    end
endfunction