function [Q, R] = qr_GS(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for j = 1:n
        coluna = A(:, j); //jº coluna de A

        //seleciona todas colunas já ortogonormalizadas
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            //gram schmidt, v_i já está normalizado
            coluna = coluna - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(coluna, 2);
        Q(:, j) = coluna/R(j, j); //normaliza a coluna ortogonal resultante
    end
endfunction

