function [lambda, x1, k, n_erro] = metodo_potencia_v1(A, x0, epsilon, M)
    k = 0;
    x0 = x0 + 1e-5;
    [valor, posicao] = max(abs(x0)); // Acha a posicão da coordenada de maior módulo
    x0 = x0 / x0(posicao); // Normaliza x0 pela coordenada de maior módulo
    x1 = A * x0; // Aproximação do autovetor dominante
    n_erro = epsilon + 1; // Inicializa o erro com um valor que entre no loop
    
    // se a primeira iteração retornar um erro menor, não executa o loop
    while k <= M && n_erro >= epsilon
        [v, p] = max(abs(x1));
        lambda = x1(p); // Aproximação do autovalor dominante
        x1 = x1 / lambda; // Normaliza x1
        n_erro = norm(x1 - x0, 'inf'); // Calcula o erro
        x0 = x1;
        x1 = A * x0;
        k = k + 1;
    end

    if k > M
        disp('O método da potência não convergiu dentro do número máximo de iterações.');
    else
        disp('O método da potência convergiu com sucesso.');
    end

    [v, p] = max(abs(x1));
    lambda = x1(p); // Autovalor dominante
    x1 = x1 / lambda; // Autovetor unitário correspondente a lambda
endfunction
