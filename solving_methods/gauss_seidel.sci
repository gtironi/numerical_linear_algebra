//////////////////////////////////////////////////////////////////////////
// Variáveis de saída:
// xk: solução do sistema Ax=b pelo método iterativo.
// norm_dif: norma da diferença entre as duas últimas iterações.
// iter: número da iteração em que o modelo atingiu a tolerância proposta.
// norm_res: norma do resíduo do sistema, ou seja, a norma de b - Axk.
//////////////////////////////////////////////////////////////////////////
function [xk, norm_dif, iter, norm_res] = gauss_seidel_inv(A, b, x0, E, M, norm_type)
    //separa A em L, D e U, com A = L + D + U
    D = diag(diag(A));
    U = triu(A - D); 
    L = tril(A - D);

    //calcula a inversa da matriz diagonal D + L, usando inv()
    inv_DL = inv(D + L);

    xk = x0;
    iter = 0;
    norm_dif = E + 1;

    //define o número maximo de iterações
    while iter <= M && norm_dif > E
        xk_old = xk; //armazena o valor anterior de xk
        xk = inv_DL * (b - U * xk_old); //atualiza xk de acordo com o método de Gauss-Seidel (inversa)

        norm_dif = norm(xk - xk_old, norm_type); // calcula a norma da diferença
        iter = iter + 1;
    end
    // calcula a norma do resíduo do sistema
    norm_res = norm(b - A*xk, norm_type);
end

//////////////////////////////////////////////////////////////////////////
// Variáveis de saída:
// xk: Solução do sistema linaer
// Obs: A matrix L de entrada deve ser triangular inferior
//////////////////////////////////////////////////////////////////////////
function x = linear_solver(L, b)
        n = length(b);
        x = zeros(n, 1);
        
        x(1) = b(1) / L(1, 1); //calcula o primeiro elemento de x diretamente
        for i = 2:n
            x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i); //define os outros elementos do vetor x
            // seleciona o valor em b, então subtrai o vetor correspondente aos elementos 
            // da linha abaixo da diagonal multiplicado pelos valores já encontrados de y
        end
    end

//////////////////////////////////////////////////////////////////////////
// Variáveis de saída:
// xk: solução do sistema Ax=b pelo método iterativo.
// norm_dif: norma da diferença entre as duas últimas iterações.
// iter: número da iteração em que o modelo atingiu a tolerância proposta.
// norm_res: norma do resíduo do sistema, ou seja, a norma de b - Axk.
//////////////////////////////////////////////////////////////////////////
function [xk, norm_dif, iter, norm_res] = gauss_seidel_sislin(A, b, x0, E, M, norm_type)
    //separa A em L, D e U, com A = L + D + U
    D = diag(diag(A));
    U = triu(A - D); 
    L = tril(A - D);

    LD = D + L;

    xk = x0;
    iter = 0;
    norm_dif = E + 1;

    while iter <= M && norm_dif > E
        xk_old = xk; //armazena o valor anterior de xk
        xk = linear_solver(LD, b -U * xk_old); //atualiza xk de acordo com o método de Gauss-Seidel
        norm_dif = norm(xk - xk_old, norm_type); // calcula a norma da diferença
        iter = iter + 1;
    end
    // calcula a norma do resíduo do sistema
    norm_res = norm(b - A*xk, norm_type);
end

