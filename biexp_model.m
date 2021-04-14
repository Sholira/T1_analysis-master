function [y_fitted, params, y_TC, resnorm] = biexp_model(dt, y, varargin)
    ylow = min(y);
    yhigh = max(y);
    contrast = (yhigh - ylow);
    
    % Initial values
    % Get the initial values if there were given to the function
    n_arg = length(varargin);
    if n_arg == 0
        x0 = [ylow, contrast/2, 1,  contrast/2, 1];
        disp('Default x0')
    else
        x0 = [varargin{:}];
    end
    
    % Parameters to find
    % x(1) = B
    % x(2) = A1
    % x(3) = t1
    % x(4) = A2
    % x(5) = t2
    % Function model
    F = @(x,dt)x(1) + x(2)*exp(-dt/x(3)) + x(4)*exp(-dt/x(5));

    % Parameters boundaries
    lb = [ylow, 1e-3, 1e-3, 1e-3, 1e-3];
    ub = [yhigh, contrast , 1e5, contrast , 1e5];

    % Fitting
    options = optimoptions('lsqcurvefit','FunctionTolerance',1e-9,'MaxIterations',1e4,'MaxFunctionEvaluations',1e4,'StepTolerance',1e-9);
    [params,resnorm,~,~,~] = lsqcurvefit(F,x0,dt,y,lb,ub,options);
    
    % Choose the slower time constant
    y_fitted = F(params,dt);
    y_TC = F(params,max([params(3) params(5)]));
    
    % Sort the parameters so the slower constant is in the fifth position
    if params(3) > params(5)
        tmp = [params(2) params(3)];
        params([2 3]) = [params(4) params(5)];
        params([4 5]) = tmp;
    end
end
   
