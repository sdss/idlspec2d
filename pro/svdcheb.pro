function svdcheb,X,M
;
;       Chebyshev polynomial basis function
;       generator for SVDFIT
;
        XX=X[0]                   ; ensure scalar XX
        sz=reverse(size(XX))      ; use size to get the type
        IF sz[n_elements(sz)-2] EQ 5 THEN $
                basis=DBLARR(M) else basis=FLTARR(M)
;
;       Calculate and return the basis functions
;
        basis[0]=1.0
        IF M ge 2 THEN basis[1]=XX
        FOR i=2,M-1 DO BEGIN
            basis[i] = 2.0 * XX * basis[i-1] - basis[i-2] 
        ENDFOR
        return,basis
end

