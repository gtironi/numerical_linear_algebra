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
X_train = dados_treinamento(:, [1,2,3,4,5,8]); //seleção de apenas algumas features
y_train = dados_treinamento(:, 11);

// Adicionar a coluna de 1s para o bias
X_train = [ones(size(X_train, 1), 1), X_train];

// Calcular os coeficientes usando minimos quadrados
weights = minimos_quadrados(X_train, y_train);

// Carregar os dados de teste
dados_teste = csvRead('data/cancer_test_2024.csv', ';');
X_test = dados_teste(:, [1,2,3,4,5,8]); //seleção de apenas algumas features
y_test = dados_teste(:, 11);

// Adicionar a coluna de 1s para o bias
X_test = [ones(size(X_test, 1), 1), X_test];

// Classificar os dados de treinamento e teste
y_train_pred = classificar(X_train, weights);
y_test_pred = classificar(X_test, weights);

// Acurácia
acuracia_treinamento = sum(y_train == y_train_pred) / length(y_train);
acuracia_teste = sum(y_test == y_test_pred) / length(y_test);

// Matriz de confusao
matriz_confusao = calcular_matriz_confusao(y_test, y_test_pred);

[precisao, recall, fpr, fnr] = calcular_metricas(matriz_confusao);

// Exibir os resultados
disp("Matriz de Confusão:");
disp(matriz_confusao);
disp("-------------------------------");
disp("Acurácia no conjunto de treinamento: " + string(acuracia_treinamento));
disp("Acurácia no conjunto de teste: " + string(acuracia_teste));
disp("-------------------------------");
disp("Precisão: " + string(precisao));
disp("Recall: " + string(recall));
disp("Probabilidade de Falso Alarme (FPR): " + string(fpr));
disp("Probabilidade de Falsa Omissão (FNR): " + string(fnr));

