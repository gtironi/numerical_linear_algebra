//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//X: as colunas são soluções do sistema linear Axi = bi, na mesma ordem de B
//////////////////////////////////////////////////////////////////////////
function [X]=Resolve_com_LU(C, B, varargin)
[n]=size(C,1);
[n_b] = size(B,2); //conta quantos valores de b estão sendo passados na matriz B

if nargin < 3 //detecta a matriz de permutação P foi passada
    P = eye(n,n); //caso não tenha sido, define P como a identidade
else
    P = varargin(1); //caso tenha, define P como a matriz passada no 3 elemento
end

L = tril(C, -1) + eye(C); //define a matriz L a partir de C, note que a diagonal é definida como 1.
U = triu(C); //define a matriz U

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