PRO defsymbols

  ;; Define symbols in true-type and vector type fonts
  ;; System variables !TSYM AND !VSYM, structures, 
  ;; are created.  
  ;;
  ;; !CSYM is also created and will contain the common symbols.
  ;; It will copy in symbols from !TSYM or !VSYM depending
  ;; on !p.font (1 or 0 for true-type, -1 for vector)
  ;; 
  ;; Thus if you change font types be sure to rerun defsymbols

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; true type font symbols (vector below)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defsysv, '!TSYM', exist=exist
  IF NOT exist THEN BEGIN 
      ;; define some plotting symbols
      ;; put !7 on first one to make times font
      ;; these are common with vector drawn fonts
      alpha = '!9a!X'
      beta = '!9b!X'
      gamma = '!9g!X'
      gamma_cap = '!9G!X'
      delta = '!9d!X'
      delta_cap = '!9D!X'
      epsilon = '!9e!X'
      zeta = '!9z!X'
      eta = '!9h!X'
      theta = '!9q!X'
      theta_cap = '!9Q!X'
      theta2 = '!9'+string(74b)+'!X'
      kappa = '!9k!X'
      lambda = '!9l!X'
      lambda_cap = '!9L!X'
      mu = '!9m!X'
      nu = '!9n!X'
      xi = '!9x!X'
      xi_cap = '!9X!X'
      pi = '!9p!X'
      pi_cap = '!9P!X'
      rho = '!9r!X'
      sigma   = '!9s!X'
      sigma_cap = '!9S!X'
      tau = '!9t!X'
      upsilon = '!9u!X'
      upsilon_cap = '!9'+string(161b)+'!X'
      phi = '!9j!X'
      phi_cap = '!9J!X'
      phi2 = '!9f!X'
      chi = '!9c!X'
      psi = '!9y!X'
      psi_cap = '!9Y!X'
      omega = '!9w!X'
      omega_cap = '!9W!X'

      angstrom = '!7'+strtrim(string(197B))+'!X'
      ;;angstrom = string(197B)

      grad = '!9'+string(209b)+'!X'
      gtequal = '!9'+string(179b)+'!X'
      ltequal  = '!9'+string(163b)+'!X'
      infinity = '!9'+string(165b)+'!X'
      notequal = '!9'+string(185b)+'!X'
      partial = '!9'+string(182b)+'!X'
      plusminus = '!9'+string(177b)+'!X'
      propto = '!9'+string(181b)+'!X'
      sqrtsym = '!9'+string(214b)+'!X'

      minus = '!9'+string(45b)+'!X'
      plus = '!9'+string(43b)+'!X'
      times = '!9'+string(180b)+'!X'

      prime = '!9'+string(162b)+'!X'
      primeprime = '!9'+string(178b)+'!X'
      degrees = '!9'+string(176b)+'!X'

      ;; not in vector sym
      upsilon_cap2 = '!9U!X'
      omega2 = '!9v!X'
      approx = '!9'+string(64b)+'!X'
      approx2 = '!9'+string(187b)+'!X'
      sum = '!9'+string(229b)+'!X'
      earth = '!9'+string(197b)+'!X'
      asterisk = '!9'+string(42b)+'!X'

      tsym1 = create_struct('alpha', alpha, $
                           'beta', beta, $
                           'gamma', gamma, $
                           'gamma_cap', gamma_cap, $
                           'delta', delta, $
                           'delta_cap', delta_cap, $
                           'epsilon', epsilon, $
                           'zeta', zeta,  $
                           'eta', eta, $
                           'theta', theta, $
                           'theta_cap', theta_cap,$
                           'theta2',theta2,$
                           'kappa', kappa, $
                           'lambda', lambda, $
                           'lambda_cap', lambda_cap, $
                           'mu', mu,$
                           'nu', nu, $
                           'xi', xi, $
                           'xi_cap', xi_cap, $
                           'pi', pi, $
                           'pi_cap', pi_cap, $
                           'rho', rho)
      tsym2 = create_struct('sigma', sigma, $
                           'sigma_cap', sigma_cap, $
                           'tau', tau, $
                           'upsilon', upsilon, $
                           'upsilon_cap', upsilon_cap,$
                            'upsilon_cap2',upsilon_cap2,$
                           'phi', phi, $
                           'phi_cap', phi_cap, $
                           'phi2', phi2,$
                           'chi', chi, $
                           'psi', psi, $
                           'psi_cap', psi_cap, $
                           'omega', omega, $
                           'omega_cap', omega_cap,$
                           'omega2', omega2, $
                           'angstrom', angstrom)

      tmath = create_struct('approx', approx, $
                           'approx2',approx2,$
                           'grad', grad,$
                           'gtequal', gtequal, $
                           'ltequal', ltequal,$
                           'infinity', infinity,$
                           'notequal', notequal,$
                           'partial', partial,$
                           'plusminus',plusminus,$
                           'propto', propto, $
                           'sqrt', sqrtsym, $
                           'sum', sum,$
                           'minus',minus,$
                           'plus',plus,$
                           'times',times)
      tother = create_struct('prime', prime, $
                             'primeprime',primeprime,$
                             'asterisk',asterisk,$
                             'degrees',degrees,$
                             'earth', earth)
      tsym = create_struct(tsym1,tsym2,tmath,tother)
      defsysv, '!TSYM', tsym
  ENDIF 

  defsysv, '!VSYM', exist=exist
  IF NOT exist THEN BEGIN 
      ;; define some plotting symbols
      ;; These are common with true type
      alpha = '!7a!X'
      beta = '!7b!X'
      gamma = '!7c!X'
      gamma_cap = '!7C!X'
      delta = '!7d!X'
      delta_cap = '!7D!X'
      epsilon = '!7e!X'
      zeta = '!7f!X'
      eta = '!7g!X'
      theta = '!7h!X'
      theta_cap = '!7H!X'
      theta2 = '!9t!X'
      kappa = '!7j!X'
      lambda = '!7k!X'
      lambda_cap = '!7K!X'
      mu = '!7l!X'
      nu = '!7m!X'
      xi = '!7n!X'
      xi_cap = '!7N!X'
      pi = '!7p!X'
      pi_cap = '!7P!X'
      rho = '!7q!X'
      sigma   = '!7r!X'
      sigma_cap = '!7R!X'
      tau = '!7s!X'
      upsilon = '!7t!X'
      upsilon_cap = '!7T!X'
      phi = '!7u!X'
      phi_cap = '!7U!X'
      phi2 = '!9P!X'
      chi = '!7v!X'
      psi = '!7w!X'
      psi_cap = '!7W!X'
      omega = '!7x!X'
      omega_cap = '!7X!X'
      
      angstrom = '!3'+strtrim(string(197b))+'!X'

      grad = '!9G!X'
      gtequal = '!9b!X'
      ltequal  = '!9l!X'
      infinity = '!9$!X'
      notequal = '!9=!X'
      partial = '!9D!X'
      plusminus = '!9+!X'
      propto = '!9?!X'
      sqrtsym = '!9S!X'
      minus = '!7-!X'
      plus  = '!7+!X'
      times = '!9X!X'

      prime = "!9'!X"
      primeprime = "!9''!X"
      degrees = '!9%!X'

      ;; not in true type
      sun = '!9n!X'
      sim = '!9A!X'
      minusplus = '!9-!X'

      vsym1 = create_struct('alpha', alpha, $
                            'beta', beta, $
                            'gamma', gamma, $
                            'gamma_cap', gamma_cap, $
                            'delta', delta, $
                            'delta_cap', delta_cap, $
                            'epsilon', epsilon, $
                            'zeta', zeta,  $
                            'eta', eta, $
                            'theta', theta, $
                            'theta_cap', theta_cap,$
                            'theta2', theta2, $
                            'kappa', kappa, $
                            'lambda', lambda, $
                            'lambda_cap', lambda_cap, $
                            'mu', mu,$
                            'nu', nu, $
                            'xi', xi, $
                            'xi_cap', xi_cap, $
                            'pi', pi, $
                            'pi_cap', pi_cap, $
                            'rho', rho)
      vsym2 = create_struct('sigma', sigma, $
                            'sigma_cap', sigma_cap, $
                            'tau', tau, $
                            'upsilon', upsilon, $
                            'upsilon_cap', upsilon_cap,$
                            'phi', phi, $
                            'phi_cap', phi_cap, $
                            'phi2', phi2,$
                            'chi', chi, $
                            'psi', psi, $
                            'psi_cap', psi_cap, $
                            'omega', omega, $
                            'omega_cap', omega_cap)

      vscience = create_struct('angstrom', angstrom, $
                               'sun', sun)
      
      vmath = create_struct('grad', grad,$
                            'gtequal', gtequal, $
                            'ltequal', ltequal,$
                            'infinity', infinity,$
                            'notequal', notequal,$
                            'partial', partial,$
                            'plusminus',plusminus,$
                            'minusplus',minusplus,$
                            'propto', propto, $
                            'sqrt', sqrtsym, $
                            'minus', minus,$
                            'plus', plus, $
                            'times',times)
      tother = create_struct('prime', prime, $
                             'primeprime',primeprime,$
                             'degrees',degrees)
      vsym = create_struct(vsym1,vsym2,vscience,vmath,tother)
      defsysv, '!VSYM', vsym
  ENDIF 

  ;; now create another system variable with only the common 
  ;; symbols from vector and truetype.  This one we can always use
  ;; to get what we want

  vtag_names = tag_names(!vsym)
  ttag_names = tag_names(!tsym)
  match, vtag_names, ttag_names, mv, mt
  nmatch=n_elements(mv)

  defsysv, '!csym', exist=cexist
  IF NOT cexist THEN BEGIN 
      csym = create_struct(vtag_names[mv[0]], '')
      FOR i=1L, nmatch-1 DO csym = create_struct(csym, vtag_names[mv[i]], '') 
          
      defsysv, '!csym', csym
  ENDIF 

  ;; copy in the correct values for the font we are using
  IF (!p.font EQ 1) OR (!p.font EQ 0) THEN BEGIN
      ;; true-type fonts
      FOR i=0L, nmatch-1 DO !csym.(i) = !tsym.(mt[i])
  ENDIF ELSE BEGIN 
      ;; vector drawn fonts
      FOR i=0L, nmatch-1 DO !csym.(i) = !vsym.(mv[i])
  ENDELSE 

END 
