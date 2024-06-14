function [lambda1, x1, k, n_erro] = potencia_deslocada_inversa(A, x0, epsilon, alfa, M)
    //Limitnado o alfa
    centros = diag(A); //Centro dos discos de Gerschgorin
    raios = sum(abs(A'), 2) - abs(centros); //Raios dos discos de Gerschgorin

    limite_superior = centros + raios;
    limite_inferior = centros - raios;

    maior_limite_superior = max(limite_superior); //Maior valores real possivel
    menor_limite_inferior = min(limite_inferior); //Menor valor real possivel

    //Atualiza o alfa caso necessário
    if alfa > maior_limite_superior
        alfa = maior_limite_superior;
    elseif alfa < menor_limite_inferior
        alfa = menor_limite_inferior;
    end

    if norm(x0, 2) == 0
        x0 = x0 + epsilon;
    end

    k = 0;
    x0 = x0 / norm(x0, 2); // Normaliza x0 pela norma 2
    n_erro = epsilon + 1; // Inicializa o erro com um valor que entre no loop
    // Decomposição LU de A - alfa * I
    [L, U, P] = lu(A - alfa * eye(A));

    while k <= M && n_erro >= epsilon
        x1 = resolve_com_LU(L, U, x0, P); // Resolve o sistema com a função Resolve_com_LU
        x1 = x1 / norm(x1, 2); // Normaliza x1
        lambda = x1' * A * x1; // Quociente de Rayleigh; x1 é unitário
        if x1' * x0 < 0
            x1 = -x1; // Mantém x1 com o mesmo sentido de x0
        end
        n_erro = norm(x1 - x0, 2); // Calcula o erro
        x0 = x1;
        k = k + 1;
    end

    if k > M
        disp('O método da potência deslocada com iteração inversa não convergiu dentro do número máximo de iterações.');
    else
        disp('O método da potência deslocada com iteração inversa convergiu com sucesso.');
    end

    lambda1 = x1' * A * x1; // Autovalor de A mais próximo de alfa
endfunction




//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//X: as colunas são soluções do sistema linear Axi = bi, na mesma ordem de B
//////////////////////////////////////////////////////////////////////////
function [X]=resolve_com_LU(L, U, B, varargin)
    [n]=size(L,1);
    [n_b] = size(B,2); //conta quantos valores de b estão sendo passados na matriz B
    
    L = L + 0.0000001;
    U = U + 0.0000001;

    if nargin < 3 //detecta a matriz de permutação P foi passada
        P = eye(n,n); //caso não tenha sido, define P como a identidade
    else
        P = varargin(1); //caso tenha, define P como a matriz passada no 3 elemento
    end
    
    // aplica a matriz de permutação em B
    // já que PLUX = B, levando a LUX = PB
    // B antes e depois da permutação são chamados de B, mas são coisas diferentes.
    B = P * B;
    
    for j=1:n_b //percorre todos os vetores b em B
        y=zeros(n,1); //zera o vetor y a cada iteração
        b = B(1:n,j); //escolhe o vetor b em B, de acordo com j
        y(1)=b(1)/L(1,1); //define o primeiro elemento do vetor y
        for i=2:n 
            y(i)=(b(i)-L(i,1:i-1)*y(1:i-1)); //define os outros elementos do vetor y
            // seleciona o valor em b, então subtrai o vetor correspondente aos elementos 
            // da linha abaixo da diagonal multiplicado pelos valores já encontrados de y
        end
        B(1:n,j) = y; //substitui o valor de y na propria matriz B para economizar memória
    end
    Y = B; //define a nova matriz como Y, contendo todos os vetores y encontrados

    for j=1:n_b //percorre todos os vetores y em Y
        x=zeros(n,1); //zera o vetor x a cada iteração
        y = Y(1:n,j); //escolhe o vetor y em Y, de acordo com j
        x(n)=y(n)/U(n,n); //define o último elemento do vetor x
        for i=n-1:-1:1
            x(i)=(y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i); //define os outros elementos do vetor x
            // seleciona o valor em y, então subtrai o vetor correspondente aos elementos 
            // da linha acima da diagonal multiplicado pelos valores já encontrados de x
            // então divide pelo elemento da diagonal naquela linha
        end
        Y(1:n,j) = x; //substitui o valor de x na propria matriz Y para economizar memória
    end
    X = Y; //define a nova matriz como X, contendo todos os vetores x encontrados
    endfunction