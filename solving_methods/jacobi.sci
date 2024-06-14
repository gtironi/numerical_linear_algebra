//////////////////////////////////////////////////////////////////////////
// Variáveis de saída:
// xk: solução do sistema Ax=b pelo método iterativo.
// norm_dif: norma da diferença entre as duas últimas iterações.
// iter: número da iteração em que o modelo atingiu a tolerância proposta.
// norm_res: norma do resíduo do sistema, ou seja, a norma de b - Axk.
//////////////////////////////////////////////////////////////////////////
function [xk, norm_dif, iter, norm_res] = jacobi(A, b, x0, E, M, norm_type)
    //separa A em L, D e U, com A = L + D + U
    D = diag(diag(A));
    U = triu(A - D); 
    L = tril(A - D);

    //calcula a inversa da matriz diagonal D
    //como D é diagonal, sua inversa é simplesmente a inversa de cada elemento na diagonal (inv(A)nn = 1/(A)nn)
    inv_D = diag(diag(1./D));
    xk = x0;
    iter = 0;
    norm_dif = E + 1;

    //define o número maximo de iterações
    while iter <= M && norm_dif > E
        xk_old = xk; //armazena o valor anterior de xk
        xk = inv_D * (b - (L + U) * xk_old); //atualiza xk de acordo com o método de Jacobi

        norm_dif = norm(xk - xk_old, norm_type); // calcula a norma da diferença
        iter = iter + 1;
    end
    // calcula a norma do resíduo do sistema
    norm_res = norm(b - A*xk, norm_type);
end

