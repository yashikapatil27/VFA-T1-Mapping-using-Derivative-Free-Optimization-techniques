% Start fresh
clear all; close all;

% Objective function
f = @(x) sin(pi*x(1)).*sin(pi*x(2));

% Gradient of Objective function
gradf = @(x)[pi.*cos(pi.*x(1)).*sin(pi.*x(2));pi.*cos(pi.*x(2)).*sin(pi.*x(1))];

% Surf and Contour plot of the objective function
figure('Position',[30 100 1200 500])
N = 60;
x = linspace(-1,1,N);
[X,Y] = meshgrid(x,x);
fplot = @(x,y) sin(pi*x).*sin(pi*y);
Z = fplot(X,Y);
subplot(1,3,1)
surf(X,Y,Z);
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
hold on;
subplot(1,3,2)
contour(X,Y,Z,60)
hold on
xlabel('x'); ylabel('y');
subplot(1,3,3)
h1 = plot(NaN,NaN,'-b');hold on; h2 = plot(NaN, NaN,'--r');
h3 = plot(NaN,NaN,'*k','MarkerSize',10);h4 = plot(NaN,NaN, '*r','MarkerSize',10);
xlabel('alpha'); ylabel('f(x,y)');
ylim([min(Z(:))-0.2 max(Z(:))+0.2]);
title('Wolfe Conditions');

% Set initial point
xk = [-0.25;-0.5];
fk = f(xk);

fprintf('iter                   xk                          fk                      alphak\n')
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',0,xk(1),xk(2),fk, 1)

% Perform iteration
c1 = 0.25; % Armijo condition constant
c2 = 0.9;  % Curvature condition constant

k = 0;
keepIterate = true;
while keepIterate
    k = k+1;
    subplot(1,3,2)
    plot(xk(1),xk(2),'o','MarkerFaceColor','r')

    %choose pk
    gradfk = gradf(xk);
    pk = -gradfk;

    %Plot wolfe curves
    alpha = linspace(0,1,121);
    
    % Plotting Wolfe curves
    xx = bsxfun(@plus,xk,bsxfun(@times,alpha,pk));
    set(h1,'Xdata',alpha,'YData',fplot(xx(1,:),xx(2,:)));
    set(h2,'Xdata',alpha,'YData',fplot(xk(1),xk(2)) + (c1*gradfk.'*pk)*alpha);

    % Backtracking linesearch using Wolfe Conditions
    alpha_star = line_search_algorithm(f, gradf, xk, pk, c1, c2);

    subplot(1,3,2)
    quiver(xk(1),xk(2),alpha_star*pk(1),alpha_star*pk(2),1,'k','Linewidth',2);

    % Step xk to xk + alpha*pk
    xk = xk + alpha_star * pk;
    fktest = f(xk); % Update function value
    funcdiff = fk-fktest;
    fk = fktest;

    fprintf('%3d        %3.8e    %3.8e       %3.8e       %3.8e\n',k,xk(1),xk(2),fk,alpha_star)
    
    %Stopping criteria
    if funcdiff < 1e-6
        disp('Stopping criteria triggered');
        break;
    end

end

function alpha_star = line_search_algorithm(f, gradf, xk, pk, c1, c2)
    alpha0 = 0;
    alpha_max = 1;
    alpha1 = 0.6;
    i = 1;

    while true
        % Evaluate phi(alpha_i)
        phi_0 = f(xk);
        phi_prime_0 = gradf(xk)' * pk;
        phi_alpha = f(xk + alpha1 * pk);

        if phi_alpha > phi_0 + c1 * alpha1 * phi_prime_0 || (phi_alpha >= phi_0 && i > 1)
            fprintf("First loop.")
            alpha_star = zoom(f, gradf, xk, pk, alpha0, alpha1, c1, c2);
            break;
        end

        % Evaluate phi_prime(alpha_i)
        phi_prime_alpha = gradf(xk + alpha1 * pk)' * pk;

        if abs(phi_prime_alpha) <= -c2 * phi_prime_0
            fprintf("Second loop.")
            %alpha_star = alpha1;
            break;
        end

        if phi_prime_alpha >= 0
            fprintf("Third loop.")
            alpha_star = zoom(f, gradf, xk, pk, alpha1, alpha0, c1, c2);
            break;
        end

        alpha0 = alpha1;
        alpha1 = min(2 * alpha1, alpha_max);
        fprintf("Continue iteration.")
        i = i + 1;
    end
end

function alpha_star = zoom(f, gradf, xk, pk, alpha_lo, alpha_hi, c1, c2)
    while true
        alpha_j = interpolate(alpha_lo, alpha_hi); % Interpolate to find a trial step length alpha_j between alpha_lo and alpha_hi
        phi_0 = f(xk);
        phi_prime_0 = gradf(xk)' * pk;
        phi_alpha_j = f(xk + alpha_j * pk);

        if phi_alpha_j > phi_0 + c1 * alpha_j * phi_prime_0 || phi_alpha_j >= f(xk + alpha_lo * pk)
            alpha_hi = alpha_j;
        else
            phi_prime_alpha_j = gradf(xk + alpha_j * pk)' * pk;

            if abs(phi_prime_alpha_j) <= -c2 * phi_prime_0
                alpha_star = alpha_j;
                break;
            end

            if phi_prime_alpha_j * (alpha_hi - alpha_lo) >= 0
                alpha_hi = alpha_lo;
            end

            alpha_lo = alpha_j;
        end
    end
end

function alpha_j = interpolate(alpha_lo, alpha_hi)
    %(quadratic, cubic, or bisection)
    alpha_j = (alpha_lo + alpha_hi) / 2; % For simplicity, use bisection
end
