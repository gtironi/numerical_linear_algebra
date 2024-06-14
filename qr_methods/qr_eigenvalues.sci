function S = espectro(A, tol)
    n = size(A, 1);
    S = diag(A);
    diff = tol + 1; //garante a entrada no loop

    //itera até a tolerância ser alcançada
    while diff > tol
        [Q, R] = qr(A); //acha a decomposição qr
        A = R * Q;
        S_new = diag(A);
        diff = norm(S_new - S, 'inf');
        S = S_new;
    end
endfunction


