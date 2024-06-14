exec('gaussian_elimination.sci');

//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//alpha: solução do sistema A_tAx=A_tb (assumimos que tal solução existe).
//////////////////////////////////////////////////////////////////////////
function alpha = minimos_quadrados(X, b)
        
    // Calcular (X^T * X)
    X_transp_X = X' * X;
    
    // Calcular (X^T * y)
    X_transp_b = X' * b;
    
    alpha = gaussian_elimination(X_transp_X, X_transp_b);
    
    endfunction

