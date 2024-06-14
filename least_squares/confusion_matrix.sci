//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//matriz_confusao: [tp, fp;
//                  fn, tn]
//onde:
//tp = Verdadeiros Positivos (True Positives)
//fp = Falsos Positivos (False Positives)
//fn = Falsos Negativos (False Negatives)
//tn = Verdadeiros Negativos (True Negatives)
//////////////////////////////////////////////////////////////////////////
function matriz_confusao = calcular_matriz_confusao(y_true, y_pred)
    tp = sum((y_true == 1) & (y_pred == 1));
    tn = sum((y_true == -1) & (y_pred == -1));
    fp = sum((y_true == -1) & (y_pred == 1));
    fn = sum((y_true == 1) & (y_pred == -1));
    
    matriz_confusao = [tp, fp; fn, tn];
endfunction

//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//precisao: Precision (precisao)
//recall: Recall (sensibilidade)
//fpr: False positive rate (probabilidade de falso alarme)
//fnr: False negative rate (probabilidade de falsa omissão de alarme)
//////////////////////////////////////////////////////////////////////////
function [precisao, recall, fpr, fnr] = calcular_metricas(matriz_confusao)
    tp = matriz_confusao(1, 1);
    fp = matriz_confusao(1, 2);
    fn = matriz_confusao(2, 1);
    tn = matriz_confusao(2, 2);
    
    precisao = tp / (tp + fp);
    recall = tp / (tp + fn);
    fpr = fp / (fp + tn);
    fnr = fn / (fn + tp);
endfunction



