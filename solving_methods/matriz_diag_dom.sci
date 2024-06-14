//////////////////////////////////////////////////////////////////////////
// Variáveis de saída:
// b: Vetor b para x sendo um vetor unitário, com os três primeiros elementos igual a 2 [2 ; 2; 2; 1; ... ; 1]
// A: Matriz quadrada nxn com valores aleatórios entre 0 e 10, com a diagonal estritamente dominante.
// Obs: a garantia da dominancia se da pela soma absoluta de todos
// os elementos da linha na diagonal daquela linha
//////////////////////////////////////////////////////////////////////////
function [A, b] = matriz_diag_dominante(n)
    // gera uma matriz A com valores aleatórios entre 0 e 10
    A = rand(n, n) * 10;

    // calcula a soma dos valores absolutos de cada linha
    soma_linhas = sum(abs(A), 2);

    // adiciona uma diagonal estritamente dominante à matriz A
    A = A + diag(soma_linhas);

    // gerar vetor x e, a partir dele, o vetor b
    x = ones(n, 1);
    x(1:3) = 2;
    b = A*x
endfunction