function [En] = function_pqrec2(d, a, V0, N, Lb, n_modes)


    m_e = 9.1091e-31;
    m_eff = .067*m_e;
    e = 1.602176565e-19;
    h_bar = 6.626e-34/(2*pi);   % m_eff dans l' AsGa (1.08

    E0 = h_bar^2*pi^2/(2*m_eff*a^2)/e*1e3;  % Mode 1 pour le puit infini en meV
    v0 = V0/E0;   % V0 en meV puis normalisation

  
    pqrec = @(x) v0 * (abs(x)>(a+d/2)/a | abs(x) <= d/2/a);

    delta=Lb/N;xb=-Lb/2+Lb/N*(0:N);

    
    vn = pqrec(xb);
    %plot(xb, vn)


    am2=-1/12/delta^2;ap2=am2;am1=4/3/delta^2;ap1=am1;A0=-5/2/delta^2;

    options.disp = 0;
    ee = ones(N+1,1);
    Lap = spdiags([-ee/12 4*ee/3 -5*ee/2 4*ee/3 -ee/12],[-2 -1 0 1 2],N+1,N+1);    %spdiags = sp pour sparse (matrice creuse)
    
    A = -1 /(pi^2 * delta^2) * Lap + spdiags(vn.',0,N+1,N+1);

    [psi, En] = eigs(A, n_modes, 'sm', options);

    En = E0 * sort(diag(En));
endfunction


