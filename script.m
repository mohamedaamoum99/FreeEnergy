function mainFunc()

    % Reactants
    intput = { 'C (CI)' ; 'S' ; 'H2O' ; 'H2' ; 'N2' ; 'O2' };
    % Initial Mol Flow ( kmol/hr )
    iniMF = [ 0.085754725 0.070024953 0.139658757 0.08383186 1.074991586 0.0618914 ];


    % Products
    output = { 'N2' ; 'O2' ; 'C (PureSolid)' ; 'CO' ; 'CO2' ; 'H2' ; 'CH4' ; 'H2O' ; 'H2S' ; 'SO2' ; 'SO3' ; 'NH3' ; 'NO' ; 'NO2' };
    % Initial Solution
    iniSL = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 ];

    % Constraints 
    Aeq = [
           0   0    1    1   1  0  1  0   0  0   0  0   0  0               % C balance
           0   0    0    0   0  1  1  1   1  0   0  2   0  0               % H balance
           0   1    0    1   2  0  0  1   0  1   1  0   1  1               % O balance
           1   0    0    0   0  0  0  0   0  0   0  2   1  1               % N balance
           0   0    0    0   0  0  0  0   1  1   1  0   0  0               % S balance
          ];
    beq = [
           4 * iniMF(1)                       % C balance
           2 * iniMF(3) + 6 * iniMF(4)        % H balance
           12 * iniMF(3) + iniMF(6)           % O balance
           3 * iniMF(5)                       % N balance
           3 * iniMF(2)                       % S balance
          ];

    % Solutions Must Be Great Than or equal 0
    LB = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];     

    % Resolution
    [x] = fmincon( @func , iniSL , [] , [] , Aeq , beq , LB , [] , [] , optimset('Algorithm','sqp', 'TolFun', 1e-6) );

    % Solutions
    for i=1:14
        fprintf('%s : %.4g kmol/hr \n', output{i}, x(i));
    end

end


% Calculate objectif function
function funObj = func(x)

    % Enthalpie libre de formation de chaque constituant dans les cdts de réf en kJ·mol-1
    Gf = [ 0 0 0 -137.15 -394.37 0 -50.49 -228.572 -33.44 -300.12 -370.95 -16.4 86.57 51.328 ];

    
    % Tempurature (K)
    T = 700 + 273.15;

    % Constante des gaz parfaits (KJ·mol-1.K-1)
    R = 8.314e-3;

    % Pression (bar)
    P = 1.01325;

    sum = 0;
    for v=1:14
        sum = sum + x(v);
    end    

    funObj = 0;
    for v=1:14

        % Activity of gaz || a(i) = coeFug(i)*Frac(i)*(P/P°) || , coeff of fugacity supposed = 1, and total pressor = 1
        av = x(v)/sum ;
        
        % Activity of liquid || a(i) = coeAct(i)*Frac(i) || , No liquid in this model
        
        % Activity of solid || a(i) = 1 ||
        if(v == 3) 
            av = 1 ;
        end



        funObj = funObj + x(v).*(Gf(v) + (R*T)*log(av));

        
    end

end