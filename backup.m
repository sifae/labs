%% Task 12
X = @(x) linspace(eps,x,100);
Y = linspace(eps,10,100);
%
f = @(x) sin(x)./x;
intTrapz = @(x) trapz(x, f(x));
intRect  = @(x) rectangles(x, f(x));
intSimps = @(x) simpson(x, f(x));
%
Y1 = arrayfun(X,Y,'UniformOutput',false);
Z = zeros(1,100);
for i = 1:length(Y1)
  
  Z(i) = intTrapz(cell2mat(Y1(i)));
end
plot(Y,Z);

function fourierApprox(f,a,b,n)
  filename = 'II-5.gif';
  x = linspace(a,b,100);
  partSum = zeros(1,length(x));
  L = (b-a)/2;

  for i = 0:n
    basis = getFunc(i,meth);
    fi = trapz(f(x).*basis(x));
    fi = fi ./ (norm(basis(x))^2);
    partSum = partSum + fi.*basis(x);
    
    plot(x, f(x), x, partSum, 'LineWidth', 2)
    xlabel('x')
    drawnow
    
    gif(i, filename)
  end
end