exec('minimos_quadrados.sci');

P = [100, 101, 112, 122, 124, 122, 143, 152, 151, 126, 155, 159, 153, 177, 184, 169, 189, 225, 227, 223, 218, 231, 179, 240];
L = [100, 105, 110, 117, 122, 121, 125, 134, 140, 123, 143, 147, 148, 155, 156, 152, 156, 183, 198, 201, 196, 194, 146, 161];
K = [100, 107, 114, 122, 131, 138, 149, 163, 176, 185, 198, 208, 216, 226, 236, 244, 266, 298, 335, 366, 387, 407, 417, 431];

lnP = log(P); //Y = log(P)
lnL = log(L); //X_1 = log(L)
lnK = log(K); //X_2 = log(K)

// Y = log(b) + alpha * X_1 + (1 - alpha) * X_2

// A = [1, X_1, X_2]
A = [ones(lnP') lnL' lnK'];
y = lnP';

// Agora resolvemos A^T.A.(alphas) = A^T.y
alphas = minimos_quadrados(A, y)

// Extração dos parâmetros
ln_b = alphas(1);
alpha = alphas(2);
beta = alphas(3);

b = exp(ln_b);

disp("Os parametros aprendidos foram: b = " + string(b) + " e alpha = " + string(alpha));

// Função de produção estimada
disp("A função de produção estimada é P = " + string(b) + " * L ^ " + string(alpha) + " * K ^ " + string(1 - alpha));

// Calculo da produção nos anos de 1910 e 1920
L_1910 = 147;
K_1910 = 208;
L_1920 = 194;
K_1920 = 407;

P_1910_est = b * L_1910 ^ alpha * K_1910 ^ (1 - alpha);
P_1920_est = b * L_1920 ^ alpha * K_1920 ^ (1 - alpha);

disp("Produção estimada para 1910: " + string(P_1910_est) + " (Produção real: 159)");
disp("Produção estimada para 1920: " + string(P_1920_est) + " (Produção real: 231)");