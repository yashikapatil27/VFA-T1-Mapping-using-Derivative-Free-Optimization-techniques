% Start fresh
clear all; close all;

% Objective function
f = @(x) sin(pi*x(1)).*sin(pi*x(2));

% Gradient of Objective function
gradf = @(x)[pi.*cos(pi.*x(1)).*sin(pi.*x(2));pi.*cos(pi.*x(2)).*sin(pi.*x(1))];

% Surf and Contour plot of the objection function
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
h1 = plot(NaN,NaN,'-b'); hold on; h2 = plot(NaN,NaN,'--r');
h3 = plot(NaN,NaN,'*k','MarkerSize',10); h4 = plot(NaN,NaN,'*r','MarkerSize',10);
xlabel('alpha'); ylabel('f(x,y)');
ylim([min(Z(:))-0.2 max(Z(:))+0.2]);
title('Wolfe Conditions');

% Set initial point
xk = [-0.05;-0.18];
fk = f(xk);

fprintf('iter                   xk                          fk                      alphak\n')
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',0,xk(1),xk(2),fk, 1)

% Perform iteration
keepIterate = true;
k = 0;
while keepIterate
    k = k + 1;
    subplot(1,3,2)
    plot(xk(1),xk(2),'o','MarkerFaceColor','r')

    % Choose pk
    gradfk = gradf(xk); % or you can normalize this
    pk = -gradfk;
    
    
    % Choose alphak using Wolfe conditions
    alpha = 0.1; % Starting alpha
    c1 = 1e-4; % Armijo condition parameter
    c2 = 0.9; % Curvature condition parameter
    rho = 0.5; % Step size reduction factor

    xx = bsxfun(@plus,xk,bsxfun(@times,alpha,pk));
    set(h1,'Xdata',alpha,'YData',fplot(xx(1,:),xx(2,:)));
    set(h2,'Xdata',alpha,'YData',fplot(xk(1),xk(2)) + (c1*gradfk.'*pk)*alpha);

    while true
        if f(xk + alpha*pk) > f(xk) + c1*alpha*gradfk.'*pk || ...
                (alpha > 1 && f(xk + alpha*pk) >= f(xk))
            % Armijo condition or weak Wolfe condition not met
            alpha = zoom(alpha, pk, xk, f, gradf, c1, c2);
            fprintf("First loop.")
        elseif abs(gradf(xk + alpha*pk).' * pk) <= -c2 * gradfk.' * pk
            % Strong Wolfe condition met 
            fprintf("Second loop.")
            break;
        elseif gradf(xk + alpha*pk).' * pk >= 0
                % Strong Wolfe condition not met
                fprintf("Third loop.")
                alpha = zoom(alpha, pk, xk, f, gradf, c1, c2);

                set(h3,'XData',alpha,'YData',f(xk + alpha*pk));
                set(h4,'XData',alpha,'YData',f(xk) + (c1*gradfk.'*pk)*alpha);
                drawnow;
                pause(0.01)
         else
                % Continue iteration
                fprintf("Fourth loop.")
                break;
         end
        
    end

    subplot(1,3,2)
    quiver(xk(1),xk(2),alpha*pk(1),alpha*pk(2),1,'k','Linewidth',2);

    % Step xk to xk + alpha*pk
    xk = xk+alpha*pk;
    funcdiff = f(xk)-fk;
    fk = f(xk); % Update function value
    
    fprintf('%3d        %3.8e    %3.8e       %3.8e       %3.8e\n',k,xk(1),xk(2),fk,alpha)
    
    % Difference in function value stop criteria
    % if funcdiff <  1e-8
    %     keepIterate = false;
    %     disp('Difference in function value is below threshold.');
    %     break;
    % end
    
    if fk == -1
        keepIterate = false;
    end

end

function alpha = zoom(alpha_init, pk, xk, f, gradf, c1, c2)
    alpha_lo = 0;
    alpha_hi = alpha_init;
    while true
        % Interpolate to find a trial step length
        alpha_j = 0.5 * (alpha_lo + alpha_hi);
        
        % Evaluate objective function at the trial step length
        fj = f(xk + alpha_j * pk);
        fk = f(xk);
        if fj > fk + c1 * alpha_j * gradf(xk).' * pk || fj >= f(xk + alpha_lo * pk)
            alpha_hi = alpha_j;
        else
            gradfj = gradf(xk + alpha_j * pk).' * pk;
            if abs(gradfj) <= -c2 * gradf(xk).' * pk
                alpha = alpha_j;
                break;
            elseif gradfj * (alpha_hi - alpha_lo) >= 0
                alpha_hi = alpha_lo;
            end
            alpha_lo = alpha_j;
        end
    end
end
