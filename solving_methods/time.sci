exec('gauss_seidel.sci');
exec('matriz_diag_dom.sci');

function [tempo_inv, tempo_sislin, diferenca] = tempo_execucao(tamanho_matriz)
    // gera matriz A com diagonal estritamente dominante e vetor b compatível
    [A, b] = matriz_diag_dominante(tamanho_matriz);

    // vetor inicial x0 de zeros
    x0 = zeros(tamanho_matriz, 1);

    // parâmetros para o método de Gauss-Seidel
    E = 1e-2;
    M = 500;
    norm_type = 2;

    // mede o tempo de execução do método de Gauss-Seidel com inv()
    tic();
    [xk, norm_dif, iter, norm_res] = gauss_seidel_inv(A, b, x0, E, M, norm_type);
    tempo_inv = toc();

    // mede o tempo de execução do método de Gauss-Seidel com a resolução do sistema linear
    tic();
    [xk, norm_dif, iter, norm_res] = gauss_seidel_sislin(A, b, x0, E, M, norm_type);
    tempo_sislin = toc();

    diferenca = abs(tempo_inv - tempo_sislin)


    disp('Para n igual a ' + string(tamanho_matriz));
    disp('Tempo do método de Gauss-Seidel com inv(): ' + string(tempo_inv) + ' segundos')
    disp('Tempo do método de Gauss-Seidel com sistema linear: ' + string(tempo_sislin) + ' segundos')
    disp('Diferença de tempo: ' + string(diferenca) + ' segundos')
endfunction
