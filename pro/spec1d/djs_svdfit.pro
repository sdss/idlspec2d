; This is a modified version of SVDFIT from IDL 5.3, since the version
; in IDL 5.4.1 crashes when the matrix is singular!!
; Presently, this is called by FIND_NMINIMA().
; $Id$

FUNCTION DJS_SVDFIT,X,Y,M, YFIT = yfit, WEIGHTS = weights, CHISQ = chisq, $
   SINGULAR = sing, VARIANCE = variance, COVAR = covar, SIGMA=sigma, $
   Function_name = function_name, A=A,LEGENDRE=legendre,DOUBLE=double
;+
; NAME:
;   SVDFIT
;
; PURPOSE:
;   Perform a general least squares fit with optional error estimates.
;
;   This version uses the Numerical Recipies (2nd Edition) function
;   SVDFIT.  A user-supplied function or a built-in polynomial or
;   legendre polynomial is fit to the data.
;
; CATEGORY:
;   Curve fitting.
;
; CALLING SEQUENCE:
;   Result = SVDFIT(X, Y, [M])
;
; INPUTS:
;   X:   A vector representing the independent variable.
;
;   Y:   Dependent variable vector.  This vector must be same length 
;      as X.
;
; OPTIONAL INPUTS:
;
;   M:   The number of coefficients in the fitting function.  For 
;        polynomials, M is equal to the degree of the polynomial + 1.
;        If not specified and the keyword A is set, then
;        M = N_ELEMENTS(A).
;
; INPUT KEYWORDS:
;
;   A:  The inital estimates of the desired coefficients. If M
;       is specified, then A must be a vector of M elements. 
;       If A is specified, then the input M can be omitted and
;       M=N_ELEMENTS(A). If not specified, the initial value
;       of each coefficient is taken to be 1.0. If both M and A are
;       specified, them must agree as to the number of paramaters.
;
;   DOUBLE:   Set this keyword to force double precision computations. This
;       is helpful in reducing roundoff errors and improves the chances
;       of function convergence.
;
;   WEIGHTS:   A vector of weights for Y[i].  This vector must be the same
;       length as X and Y.  If this parameter is ommitted, 1's
;       (No weighting) are assumed.  
;	The error for each term is weighted by Weight[i] when computing the
; 	fit.  Gaussian or instrumental uncertianties should be weighted as
;	Weight = 1/Sigma where Sigma is the measurement error or standard
; 	deviations of Y. For Poisson or statistical weighting use
; 	Weight=1/sqrt(Y), since Sigma=sqrt(Y).
;
; FUNCTION_NAME:
;      A string that contains the name of an optional user-supplied 
;      basis function with M coefficients. If omitted, polynomials
;      are used.
;
;      The function is called: R=SVDFUNCT(X,M)
;
;      where X and M are  scalar values, and the function value is an 
;      M element vector evaluated at X with the M basis functions.  
;               M is the degree of the polynomial +1 if the basis functions are
;               polynomials.  For example, see the function SVDFUNCT or SVDLEG,
;               in the IDL User Library:
;
;      For more examples, see Numerical Recipes in C, second Edition,
;               page 676-681.
;
;      The basis function for polynomials, is R[j] = x)^j.
;
;           The function must be able to return R as a FLOAT vector or
;               a DOUBLE vector depending on the input type of X.
;      
;     LEGENDRE: Set this keyword to use the IDL function SVDLEG in the lib
;               directory to fit the data to an M element legendre polynomial.
;               This keyword overrides the FUNCTION_NAME keyword.
;
; OUTPUTS:
;   SVDFIT returns a vector of the M coefficients fitted to the
;   supplied function.
;
; OPTIONAL OUTPUT PARAMETERS:
;
;   CHISQ:   Sum of squared errors multiplied by weights if weights
;      are specified.
;
;   COVAR:   Covariance matrix of the coefficients.
;
;    VARIANCE:   Sigma squared in estimate of each coeff(M).
;               That is sqrt(VARIANCE) equals the 1 sigma deviations
;               of the returned coefficients.
;
;      SIGMA:   The 1-sigma error estimates of the returned parameters,
;               SIGMA=SQRT(VARIANCE).
;
;    SINGULAR:   The number of singular values returned.  This value should
;      be 0.  If not, the basis functions do not accurately
;      characterize the data.
;
;   YFIT:   Vector of calculated Y's.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   None.
;
; MODIFICATION HISTORY:
;   Adapted from SVDFIT, from the book Numerical Recipes, Press,
;   et. al., Page 518.
;   minor error corrected April, 1992 (J.Murthy)
;
;   Completely rewritten to use the actual Numerical Recipes routines
;   of the 2nd Edition (V.2.06). Added the DOUBLE, SIGMA, A, and
;   LEGENDRE keywords. Also changed Weight to Weights to match the
;   other fitting routines.
;
;       EXAMPLE:
;-
   ON_ERROR,2      ;RETURN TO CALLER IF ERROR
;
;       Ensure that the proper number/type of parameters are used
;
        np=N_PARAMS()
        HAVE_IN_A=KEYWORD_SET(A)
        IF HAVE_IN_A THEN IN_A=A
        IF NOT KEYWORD_SET(LEGENDRE) THEN use_legendre=0 ELSE use_legendre=1

;
;       Either two or three parameters must be supplied
;
        IF np NE 2 and np NE 3 THEN $
     MESSAGE,'Incorrect number of arguments for SVDFIT,x,y,[m]'

;
;       If 2 paramaters are supplied, then A must be supplied
;
        IF np EQ 2 THEN $
          IF N_ELEMENTS(A) eq 0 THEN $
        MESSAGE,'The keyword A must be set to call SVDFIT as SVDFIT,x,y' $
          ELSE M=N_ELEMENTS(A)
;
;
   IF np EQ 3 THEN BEGIN
     IF (SIZE(M))[0] NE 0 THEN $
           MESSAGE,'The input M must be a scalar.' 
     IF not KEYWORD_SET(A) THEN A=REPLICATE(1.0,M) ELSE $
        IF N_ELEMENTS(A) NE M or (SIZE(A))[0] NE 1 THEN $ 
           MESSAGE,'The keyword A must be an M element vector.' 
        ENDIF

   IF (SIZE(X))[0] NE 1 THEN $
           MESSAGE,'The input X must be a vector.' 

   IF (SIZE(Y))[0] NE 1 THEN $
           MESSAGE,'The input Y must be a vector.' 
    
;
;       Evaluate/Set/Check a few essential paramaters
;

   THRESH = 1.0E-9      ;Threshold used in editing singular values

   ndata = N_ELEMENTS(x)    ;SIZE of X

   IF ndata NE N_ELEMENTS(y) THEN  $ 
     message, 'Error: X and Y must have same # of elements.'

   IF N_ELEMENTS(weights) NE 0 THEN BEGIN
      IF N_ELEMENTS(weights) NE ndata THEN $
        MESSAGE, 'Error: Weights must have the number of elements as X and Y.'
     sig = 1.0/weights          ;Apply weights
   ENDIF ELSE sig = REPLICATE(1.0, ndata)

   IF n_elements(FUNCTION_NAME) EQ 0 THEN FUNCTION_NAME='svdfunct' 
   IF FUNCTION_NAME EQ '' THEN FUNCTION_NAME='svdfunct' 
      IF use_legendre THEN FUNCTION_NAME='svdleg'

;
;       If x or y or sig or a is double precision, set use_double to true
;
; The following line edited to agree with the IDL 5.4.1 code (DJS 11-Feb-2002)
;      IF ARG_PRESENT(double) then use_double = KEYWORD_SET(double) $
      IF (N_ELEMENTS(double) GT 0) then use_double = KEYWORD_SET(double) $
      else use_double = (reverse(size(a)))[1] eq 5 or $
        (reverse(size(x)))[1] eq 5 or $
        (reverse(size(y)))[1] eq 5 or $
        (reverse(size(sig)))[1] eq 5

      IF use_double THEN begin
          chisq=0.0D
          xx=double(x) & yy=double(y) & a=double(a)
          ssig=abs(double(sig)) > 1E-12
      ENDIF ELSE BEGIN
          chisq=0.0
          xx=float(x) & yy=float(y) & a=float(a)
          ssig=abs(float(sig)) > 1E-6
      ENDELSE

;
;       Call the actual NR routine.
;   
;       Warning, direct use of this function is not supported, as
;       the calling sequence is scheduled to change in a future release
;       of IDL.
;
        NR__SVDFIT,FUNCTION_NAME,xx,yy,ssig,ndata,A,M,COVAR,CHISQ,$
         DOUBLE=use_double
      
;
;       The variance is the diagonal elements of the covariance matrix
;
        variance=COVAR[lindgen(M)*(M+1)];
        sigma=sqrt(abs(variance))

        small=WHERE(variance LE max(variance)*thresh, cc)
        IF cc GT 0 THEN variance[small]=0

   good = WHERE(variance GT 0, ng) ;Cutoff for sing values
   sing = M - ng      ;# of singular values
   IF sing NE 0 THEN BEGIN
      message, 'Warning: ' + strcompress(sing, /REMOVE) + $
         ' singular values found.',/INFORM
;
      IF ng EQ 0 THEN return,undefined
   ENDIF
;
;     calculate YFIT 
;
         IF KEYWORD_SET(use_double) THEN BEGIN
            yfit=DBLARR(ndata)
            afunc=DBLARR(M)
         ENDIF ELSE BEGIN
            yfit=FLTARR(ndata)
            afunc=FLTARR(M)
         ENDELSE
; The following line edited to agree with the IDL 5.4.1 code (DJS 11-Feb-2002)
         FOR i=0L,ndata-1 DO $
            YFIT[i]=total(A * call_function(FUNCTION_NAME, xx[i], M))
;
;     Return the fitted parameters
;
         OUT_A=A
         IF HAVE_IN_A THEN A=IN_A
    RETURN,OUT_A
end
