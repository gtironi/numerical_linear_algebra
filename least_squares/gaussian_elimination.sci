//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//x: solução do sistema Ax=b (assumimos que tal solução existe).
//D: Seja A=LU a decomposição LU de A.
//Então C(i,j)=L(i,j) para i>j e C(i,j)=U(i,j) para j>=i.
//P: é a matriz de permutação de PA = LU
//////////////////////////////////////////////////////////////////////////
function [x, D, P]=gaussian_elimination(A, b)
[n]=size(A,1);
I = eye(n, n); //define uma matriz identidade n x n
C=[A,b,I]; //junta as matrizes A, o vetor b e a matriz identidade criada
for j=1:(n-1)
 //O pivô está na posição (j,j)

 // procura o maior valor em modulo nas linhas abaixo
 // e guarda a posição dele na variavel index
 [_, index] = max(abs(C(j:n, j))); 

 index = index - 1; //corrige o valor do index, pois o vetor procurado inclui o pivô
 C([j, j+index], :) = C([j+index, j], :); //realiza a troca das linhas
 for i=(j+1):n
//O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
 C(i,j)=C(i,j)/C(j,j);
 //Linha i  Linha i - C(i,j)*Linha j
//Somente os elementos da diagonal ou acima dela são computados
//(aqueles que compõem a matrix U)
 C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
 end
end

x=zeros(n,1);
// Calcula x, sendo Ux=C(1:n,n+1)
x(n)=C(n,n+1)/C(n,n);
for i=n-1:-1:1
 x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
end
D = C(1:n,1:n); //isola a matriz que contêm a decomposição LU
P = C(1:n,n+2:n*2+1); //isola a matriz P
endfunction