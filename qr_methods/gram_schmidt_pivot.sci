function [Q, R, P] = qr_GSP(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    P = eye(n, n);

    for j = 1:n
        //verifica a coluna de maior norma
        norms = zeros(1, n-j+1);

        for k = j:n //calcula a norma de todas as colunas
            norms(k-j+1) = norm(A(:, k));
        end
        [value, idx] = max(norms); //verifica o indice da coluna de maior norma
        idx = idx + j - 1; 

        if idx ~= j //permuta as colunas j e idx
            A(:, [j, idx]) = A(:, [idx, j]);
            P(:, [j, idx]) = P(:, [idx, j]);
        end
        
        //gram schmidt modificado
        coluna = A(:, j); //jº coluna de A

        //seleciona todas colunas já ortogonormalizadas
        for i = 1:j-1
            R(i, j) = Q(:, i)' * coluna;
            //gram schmidt, v_i já está normalizado
            coluna = coluna - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(coluna, 2);
        Q(:, j) = coluna/R(j, j); //normaliza a coluna ortogonal resultante
    end
endfunction

