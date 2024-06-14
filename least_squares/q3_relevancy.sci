exec('minimos_quadrados.sci');
exec('matriz_confusao.sci');

// Função de classificação
function y_pred = classificar(X, weights)
    y_pred = X * weights; // Determina o hiperplano
    y_pred = sign(y_pred); // Classifica baseado no sinal de h(x)
endfunction

// Carregar os dados de treinamento
dados_treinamento = csvRead('data\cancer_train_2024.csv', ';');

// Separação das colunas de features e o alvo
X_train = dados_treinamento(:, 1:10);
y_train = dados_treinamento(:, 11);

// Adicionar a coluna de 1s para o bias
X_train = [ones(size(X_train, 1), 1), X_train];

// Calcular os coeficientes usando minimos quadrados
weights = minimos_quadrados(X_train, y_train);

disp("Coeficientes:");
disp(string(weights(1)) + ": bias");
disp(weights(2:11));