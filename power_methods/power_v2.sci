function [lambda, x1, k, n_erro] = metodo_potencia_v2(A, x0, epsilon, M)
    k = 0;
    x0 = x0 + 1e-5;
    x0 = x0 / norm(x0, 2); // Normaliza x0 pela norma 2
    x1 = A * x0; // Aproximação do autovetor dominante
    n_erro = epsilon + 1; // Inicializa o erro com um valor que entre no loop
    
    while k <= M && n_erro >= epsilon
        lambda = x1' * x0; // Quociente de Rayleigh
        if lambda < 0
            x1 = -x1; // Mantém x1 com o mesmo sentido de x0
        end
        x1 = x1 / norm(x1, 2); // Normaliza x1
        n_erro = norm(x1 - x0, 2); // Calcula o erro
        x0 = x1;
        x1 = A * x0;
        k = k + 1;
    end

    if k > M
        disp('O método da potência não convergiu dentro do número máximo de iterações.');
    else
        disp('O método da potência convergiu com sucesso.');
    end

    lambda = x1' * x0; // Autovalor dominante
    x1 = x1 / norm(x1, 2); // Autovetor unitário correspondente a lambda
endfunction


