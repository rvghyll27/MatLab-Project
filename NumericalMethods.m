function NumericalMethods
    equation_str = input('Enter function with right-hand side zero (in terms of x): ', 's');
    equation = str2func(['@(x)' equation_str]);

    method_choice = input('Choose method\n 1: Graphical\n 2: Incremental\n 3: Bisection\n 4: Regula Falsi\n 5: Simple Fixed Point\n 6: Newton Rapson\n 7: Secant\n Choose: ');

    if method_choice == 1
    % Graphical Method
    fprintf('\n\nGraphical Method\n');
    
    xL = input('Enter the initial value for x_L: ');
    
    fprintf('\nRoot Table:\n');
    fprintf('xL\t\tf(xL)\n');
    
    f_xL_prev = equation(xL);
    fprintf('%f\t%f\n', xL, f_xL_prev);
    
    xL = xL + 0.2;
    
    while true
        f_xL = equation(xL);
        fprintf('%f\t%f\n', xL, f_xL);
        
        % Check if the sign of f(xL) has changed
        if f_xL * f_xL_prev < 0
            break;
        end
        
        f_xL_prev = f_xL;
        xL = xL + 0.2;
    end
    
    % Find the root within the specified range
    root = fzero(equation, [xL - 0.2, xL]);

    fprintf('\nRoot found at: %f\n', root);

    % Plotting the function with the root marked
    xu_graphical = xL + 1.2;  % Set upper bound for plotting
    x_vals_graphical = linspace(xL - 1.4, xu_graphical, 1000); % Adjust the range for better visualization
    y_vals_graphical = equation(x_vals_graphical);
    plot(x_vals_graphical, y_vals_graphical, 'LineWidth', 2);
    hold on;
    scatter(root, equation(root), 100, 'ro', 'filled'); % Mark root on the graph
    title('Graphical Method');
    xlabel('x');
    ylabel('f(x)');
    grid on;
    hold off;
    

    elseif method_choice == 2
        % Incremental Method
        fprintf('\n\nIncremental Method\n');

        xl_incremental = input('Enter the lower bound for incremental method: ');
        delta_x_incremental = input('Enter the delta x for incremental method: ');

        % Initialize variables
        ideal_error = 0.001;
        max_iterations = 1000;

        fprintf('\nIterations\t xL\t\t ?x\t\t xU\t\t\t f(xL)\t\t\t f(xU)\t\t\t Error\t\t\t f(xL)*f(xU)\n');

        for i = 1:max_iterations
            xL = xl_incremental;
            xU = xL + delta_x_incremental;
            f_xL = equation(xL);
            f_xU = equation(xU);
            relative_error = abs((xU - xl_incremental) / xl_incremental) * 100;
            product_f_xL_xU = f_xL * f_xU;

            fprintf('%d\t\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', i, xL, delta_x_incremental, xU, f_xL, f_xU, relative_error, product_f_xL_xU);

            if abs(product_f_xL_xU) < eps || relative_error <= ideal_error
                break;
            end

            if product_f_xL_xU < 0
                % Change the increment to be smaller
                delta_x_incremental = delta_x_incremental / 10.0;
            else
                % Move xU to the current xL
                xl_incremental = xU;
            end
        end

        % Check if a root is found
        if abs(f_xL) < eps
            root = xL;
        elseif abs(f_xU) < eps
            root = xU;
        else
            root = xL; % Choose latest xL as root if error <= 0.001
        end

        fprintf('\nRoot is: %.6f\n', root);

        % Plotting the function with the root marked
        x_vals_incremental = xl_incremental:0.01:xU;
        y_vals_incremental = equation(x_vals_incremental);

        figure;
        plot(x_vals_incremental, y_vals_incremental, 'LineWidth', 2); % Plot the function
        hold on;
        scatter(root, equation(root), 100, 'ro', 'filled'); % Mark root on the graph
        title('Incremental Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;

    elseif method_choice == 3
        % Bisection Method
        fprintf('\n\nBisection Method\n');

        xl_bisection = input('Enter the lower bound for bisection method: ');
        xu_bisection = input('Enter the upper bound for bisection method: ');
        tol_bisection = input('Enter the allowed error for bisection method: ');
        max_iterations = 1000;

        % Bisection algorithm
        fprintf('\nXl\t\t\tXr\t\t\tXu\t\t\tf(Xl)\t\tf(Xr)\t\tTolerance Comparison\n');

        for i = 1:max_iterations
            xl = xl_bisection;
            xr = (xu_bisection + xl_bisection) / 2;
            fa = equation(xl);
            fc = equation(xr);

            % Calculate relative error
            relative_error = abs((xr - xl_bisection) / xr) * 100;

            fprintf('%f\t%f\t%f\t%f\t%f\t', xl_bisection, xr, xu_bisection, fa, fc);

            % Check if the relative error is below the tolerance
            if fc == 0 || relative_error < tol_bisection
                fprintf('%f (Converged)\n', relative_error);
                break;
            else
                fprintf('%f (Not Converged)\n', relative_error);
            end

            if equation(xl_bisection) * fc < 0
                xu_bisection = xr;
            else
                xl_bisection = xr;
            end
        end

        fprintf('\nRoot is: %f\n', xr);

        

        % Plotting the function with the root marked
        x_vals_bisection = linspace(xl_bisection, xu_bisection, 1000);
        y_vals_bisection = equation(x_vals_bisection);

        figure;
        plot(x_vals_bisection, y_vals_bisection, 'LineWidth', 2);
        hold on;
        plot(xr, fc, 'ro', 'MarkerSize', 10);
        title('Bisection Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;

    elseif method_choice == 4
        % Regula Falsi Method
        fprintf('\nRegula Falsi Method\n');

        xl_regula = input('Enter the lower bound: ');
        xu_regula = input('Enter the upper bound: ');
        tol_regula = input('Enter the allowed error: ');
        max_iterations = 1000;

        % Initialize the first xr
        xl = xl_regula;
        xu = xu_regula;
        fa = equation(xl);
        fu = equation(xu);
        xr = (fu * xl - fa * xu) / (fu - fa);
        fc = equation(xr);

        fprintf('\nXl\t\t\tXu\t\t\tXr\t\t\tf(Xl)\t\tf(Xr)\t\tf(xr)*f(xl)\t\tTolerance Comparison\n');

        for i = 1:max_iterations
            fa = equation(xl);
            fu = equation(xu);
            xr = (fu * xl - fa * xu) / (fu - fa);
            fc = equation(xr);

            % Calculate relative error
            if i > 1  % Avoid division by zero in the first iteration
                relative_error = abs((xr - xr_prev) / xr) * 100;
            else
                relative_error = tol_regula + 1; % Ensure the loop starts
            end

            % Determine the sign of f(xr) * f(xl)
            if fc * fa < 0
                sign = '< 0';
            elseif fc * fa > 0
                sign = '> 0';
            else
                sign = '0';
            end

            fprintf('%f\t%f\t%f\t%f\t%f\t%s\t\t\t\t%f', xl, xu, xr, fa, fc, sign, relative_error);

            % Check if the relative error is below the tolerance
            if fc * fa == 0 || relative_error < tol_regula
                fprintf('(Converged)\n');
                break;
            else
                fprintf('(Not Converged)\n');
            end

            % Update bounds
            if fa * fc < 0
                xu = xr;
            else
                xl = xr;
            end

            % Store the previous xr
            xr_prev = xr;
        end

        fprintf('\nRoot is: %f\n', xr);

       

        % Plotting the function with the root marked
        x_vals_regula = linspace(xl_regula, xu_regula, 1000);
        y_vals_regula = equation(x_vals_regula);

        figure;
        plot(x_vals_regula, y_vals_regula, 'LineWidth', 2);
        hold on;
        plot(xr, fc, 'ro', 'MarkerSize', 10);
        title('Regula Falsi Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;

    elseif method_choice == 5
    simpleFixedPointMethod(equation_str);
     




    elseif method_choice == 6
        % Newton-Raphson Method
        fprintf('\nNewton-Raphson Method\n');

        xl_newton = input('Enter the initial guess: ');
        tol_newton = input('Enter the allowed error: ');
        max_iterations = 1000;

        % Newton-Raphson algorithm
        fprintf('\nIteration\t\tx0\t\t\tf(x0)\t\t\tf''(x0)\t\t\tx1\t\t\tError\n');

        x0 = xl_newton;
        for i = 1:max_iterations
            f_x0 = equation(x0);
            df_x0 = numerical_derivative(equation, x0);

            if df_x0 == 0
                fprintf('Derivative is zero. Method fails.\n');
                break;
            end

            x1 = x0 - f_x0 / df_x0;
            error = abs((x1 - x0) / x1) * 100;

            fprintf('%d\t\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n', i, x0, f_x0, df_x0, x1, error);

            if error < tol_newton
                break;
            end

            x0 = x1;
        end

        fprintf('\nRoot is: %.9f\n', x0);


        % Plotting the function with the root marked
        x_vals_newton = linspace(x0 - 5, x0 + 5, 1000);
        y_vals_newton = arrayfun(@(x) equation(x), x_vals_newton);

        figure;
        plot(x_vals_newton, y_vals_newton, 'LineWidth', 2);
        hold on;
        plot(x0, equation(x0), 'ro', 'MarkerSize', 10);
        title('Newton-Raphson Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;

    elseif method_choice == 7
        % Secant Method
        fprintf('\nSecant Method\n');

        xl_secant = input('Enter the first guess: ');
        xu_secant = input('Enter the second guess: ');
        tol_secant = input('Enter the allowed error: ');
        max_iterations = 1000;

        % Secant algorithm
        fprintf('\nIteration\t\tx0\t\t\tx1\t\t\tf(x0)\t\t\tf(x1)\t\t\tx2\t\t\tError\n');

        x0 = xl_secant;
        x1 = xu_secant;
        for i = 1:max_iterations
            f_x0 = equation(x0);
            f_x1 = equation(x1);

            if (f_x1 - f_x0) == 0
                fprintf('Denominator is zero. Method fails.\n');
                break;
            end

            x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
            error = abs((x2 - x1) / x2) * 100;

            fprintf('%d\t\t\t\t%.4f\t\t%.4f\t\t%.4f\t\t   %.4f\t\t%.4f\t\t   %.4f\n', i, x0, x1, f_x0, f_x1, x2, error);

            if error < tol_secant
                break;
            end

            x0 = x1;
            x1 = x2;
        end

        fprintf('\nRoot is: %.4f\n', x2);


        % Plotting the function with the root marked
        x_vals_secant = linspace(x0 - 5, x0 + 5, 1000);
        y_vals_secant = arrayfun(@(x) equation(x), x_vals_secant);

        figure;
        plot(x_vals_secant, y_vals_secant, 'LineWidth', 2);
        hold on;
        plot(x2, equation(x2), 'ro', 'MarkerSize', 10);
        title('Secant Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;

    else
        fprintf('Invalid method choice.\n');
    end
end

    function simpleFixedPointMethod(equation_str)
    try
        equation = str2func(['@(x)' equation_str]);
        x0 = input('Enter the initial guess: ');
        tol = input('Enter the allowed error: ');

        max_iterations = 1000;
        iterations = 0;
        error = inf;
        data = cell(max_iterations, 4);

        while error > tol && iterations < max_iterations
            iterations = iterations + 1;
            x1 = g(x0, equation);  % Fixed-point iteration: x1 = g(x0)

            if iterations == 1
                error = inf;
            else
                error = abs((x1 - str2double(data{iterations-1, 3})) / x1) * 100;
            end

            data{iterations, 1} = iterations;
            data{iterations, 2} = num2str(x0, '%.9f');
            data{iterations, 3} = num2str(x1, '%.9f');
            data{iterations, 4} = num2str(error, '%.9f');

            x0 = x1;
        end

        root = x0;
        displayResult(data, root);

        % Plotting the function with the root marked
        x_vals_fixedpoint = linspace(root - 1, root + 1, 1000);
        y_vals_fixedpoint = arrayfun(@(x) equation(x), x_vals_fixedpoint);

        figure;
        plot(x_vals_fixedpoint, y_vals_fixedpoint, 'LineWidth', 2);
        hold on;
        plot(root, equation(root), 'ro', 'MarkerSize', 10);
        title('Simple Fixed Point Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        hold off;
    catch ME
        disp(ME.message);
    end
end


function x1 = g(x0, equation)
    % Example: Let's choose g(x) = x - f(x) / f'(x), similar to the Newton-Raphson method
    fx = evaluateFunction(equation, x0);
    fDashX = numerical_derivative(equation, x0);

    % Avoid division by zero
    if fDashX == 0
        error('Derivative is zero. Division by zero error.');
    end

    x1 = x0 - (fx / fDashX);
end

function df = numerical_derivative(func, x)
    % Numerical differentiation using central difference
    h = 1e-5; % Small step size
    df = (func(x + h) - func(x - h)) / (2 * h);
end

function f = evaluateFunction(equation, x)
    try
        f = equation(x);
    catch
        error('Error evaluating the function.');
    end
end


function displayResult(data, root)
    disp('Iteration       x0              x1              Error');
    for i = 1:size(data, 1)
        if isempty(data{i, 1})
            break;
        end
        fprintf('%d        %s        %s        %s\n', data{i, :});
    end
    fprintf('\nRoot is: %.9f\n', root);
end