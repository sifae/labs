%% Task 1
h = 0.1;
a = -2;
b = 2;
xx = (a:h/10:b);
x  = (a:h:b);
f = @(x) sign(x);

compareInterp(x,xx,f)
%% Task 2
tiledlayout(2,2)
nexttile
plotInterp(@(x) exp(x),0,2,0.2,'nearest')
hold on
plotInterp(@(x) floor(x.^2),0,2,0.05,'nearest')
title('Nearest neighbor')

nexttile
plotInterp(@(x) 2*sinc(x),-1,1,pi/11,'pchip')
hold on
plotInterp(@(x) sign(x),-1,1,0.1,'pchip')
title('Pchip')

nexttile
plotInterp(@(x) diric(x,40)+1.2,-1,1,0.1,'linear')
hold on
plotInterp(@triangleWave,-1,1,0.1,'linear')
title('Linear')

nexttile 
plotInterp(@(x) sign(x)/2+1,-2,2,0.1,'spline')
hold on
plotInterp(@(x) diric(x,20)-1,-2,2,0.1,'spline')
title('Spline')
%% Task 3 (TODO)
xmin = -5;
xmax = 5;
n_better = 100;
n_lower = 10;
xx = linspace(xmin, xmax, n_better);
x = linspace(xmin, xmax, n_lower);
f = @(t) cos(t).^2;
fmax = 2;
interp_res = interp1(x, f(x), xx, 'nearest');
h = 0.1;
apriorn = fmax * (h^4);
aposteriorn = abs(f(xx) - interp_res);
plot(xx, apriorn*ones(1, n_better), 'g-', xx, aposteriorn, 'b-');
legend('Априорная погрешность', 'Полученная погрешность');
g = @(t) t;
fmax = 1;
interp_res = interp1(x, g(x), xx, 'nearest');
h = (xmax - xmin) / n_lower;
apriorn = fmax * (h^4);
aposteriorn = abs(f(xx) - interp_res);
plot(xx, apriorn*ones(1, n_better), 'g-', xx, aposteriorn, 'b-');
legend('Априорная погрешность', 'Полученная погрешность');
%% Task 4
a = -2;
b = 2;
N = 50;
fn = @(n,x) (1+x./n).^n;
f = @exp;

% Check if parpool already exists
if isempty(gcp('nocreate'))
  parpool;
end
f = parfeval(@convergenceFunc,0,fn,f,a,b,N,'uniform');
fetchOutputs(f);
cancel(f);
%% Task 5
N = 50;
f = @(x) sign(sin(x));

% Check if parpool already exists
if isempty(gcp('nocreate'))
  parpool;
end
f1 = parfeval(@fourierApprox,0, ...
              f,-2*pi,2*pi,N,'trigonometric','files/II-5-trig.gif');
f2 = parfeval(@fourierApprox,0, ...
              f,-0.9995,0.9995,N,'chebyshev','files/II-5-cheb.gif');
f3 = parfeval(@fourierApprox,0, ...
              f,-1,1,N,'gegenbauer','files/II-5-gegenb.gif');
fetchOutputs([f1, f2, f3]);
cancel([f1, f2, f3]);
%% Task 6
f = @(x) sin(x)./x;
x = linspace(-8*pi,0,1000);

[fMax, indMax] = max(f(x));
localMin = x(islocalmin(f(x)));
[dist,idx] = min(abs(localMin - x(indMax)));
% Comet boundaries
a = x(indMax);
b = localMin(idx);
% Comet domain
if a <= b
  xx = x(x>=a & x <= b);
else 
  xx = flip(x(x<=a & x >= b));
end

plot(x,f(x),x(indMax),fMax,'r*',localMin,f(localMin),'k.')
xlabel('x')
ylabel('f(x)')
hold on
comet(xx,f(xx))
%% Task 8
cla reset
N = 500;
rhoC = @(l) supportCircle(l,[1; 2],2);
rhoR = @(l) supportRhombus(l,[2; 3],1);
rhoE = @(l) supportEllipse(l,[-4; -1],[4 0;0 1]);
rhoS = @(l) supportSquare(l,[1; 2],1);

% funs = {rhoC, rhoR, rhoE, rhoS};
% demo(funs, N);
drawSet(rhoC, N, -6, 6)
%% Task 9
N = 150;

% rhoC = @(l,a,r) supportLebesgue(@(x) norm(x-a)-r,l,a);
% rhoR = @(l,a,r) supportLebesgue(@(x) sum(abs(x-a))-r,l,a);
% rhoE = @(l,p,P) supportLebesgue(@(x) dot(x-p,P^(-1)*(x-p))-1,l,a);
% rhoS = @(l,a,r) supportLebesgue(@(x) max(abs(x-a))-r,l,a);
% 
% funs = {@(l) rhoC(l,[1; 1],3), @(l) rhoR(l,[-1; -2],2), ...
%         @(l) rhoE(l,[0; 1],[9 4;3 4]) @(l) rhoS(l,[-1; 0],1)};
%       
% demo(funs, N)
rho = @(l) supportLebesgue(@(x) x(1)+(x(2)+1)^2+2,l,[0;0]);
drawSet1(rho,N,-4,4)
%% Task 10
rhoC = @(l) supportCircle(l,[1; 0],1/3);
rhoR = @(l) supportRhombus(l,[0; 0],1);
rhoS = @(l) supportSquare(l,[0; 0], 1);
rhoE = @(l) supportEllipse(l,[0; 0], [2 0; 0 1]);

tiledlayout(2,2)
nexttile
drawPolar(rhoC,300)
nexttile
drawPolar(rhoR,300)
nexttile
drawPolar(rhoS,300)
nexttile
drawPolar(rhoE,300)
%% Functions
% Task 1
function compareInterp(x,xx,f)
  y = f(x);
  yy = zeros(4,length(xx));
  yy(1,:) = interp1(x,y,xx,'nearest');
  yy(2,:) = interp1(x,y,xx,'linear');
  yy(3,:) = interp1(x,y,xx,'spline');
  yy(4,:) = interp1(x,y,xx,'pchip');
  names = ["Nearest", "Linear", "Spline", "Pchip"];
  
  tiledlayout(2,2)
  for i = (1:4)
    nexttile
    plot(xx, f(xx), xx, yy(i,:), 'LineWidth', 2)
    xlabel('x')
    ylabel('f(x)')
    legend('Origin','Interpolation')
    title(names(i))
  end
end

% Task 2
function plotInterp(f,a,b,h,method)
  xx = (a:h/10:b);
  x  = (a:h:b);
  yy = interp1(x,f(x),xx,method);
  %Dummy points to use legend without color bindings
  plot([0 0], [0 0],'k', [0 0], [0 0],'k--'); 
  hold on
  plot(xx, yy, xx, f(xx), '--', 'LineWidth', 2)
  legend('Interpolation','Original');
end
function func = triangleWave(x)
  tmp = x - floor(x+1/2);
  func = 2.*abs(2.*tmp)-1;
end

% Task 3
function func = linearInterpError(f,xx)
  a = xx(1);
  M1 = @(x) max(diff(f(xx(xx<=x))));
  M2 = @(x) max(diff(f(xx)));
  func = @(x) (x-a).*(M1(x) + M2(x));
end

% Task 4
function convergenceFunc(fn,f,a,b,n,convType)
  filename = 'files/II-4.gif';
  x = linspace(a,b,1000);
  
  metric = getMetric(convType);
  for i = 1:n
    %Solve problem with quality and parallel processing
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(x, f(x), x, fn(i,x), 'LineWidth', 2)
    fmt = sprintf('$ %s(f_n,f) = %f $', '\rho',metric(fn(i,x),f(x)));
    title(fmt,'Interpreter','latex','FontSize',15)
    xlabel('x')
    drawnow
    
    gif(i,filename)
  end
end
function metric = getMetric(convType)
  switch convType
    case "uniform"
      metric = @(f1,f2) max(abs(f1-f2));
    case "quadratic"
      metric = @(f1,f2) sqrt(trapz((f1-f2).^2));
    case "pointwise"
      metric = @(f1,f2) missing;
    otherwise
      fmt = ['Error. \nconvType must be: ', ...
            '"pointwise", "uniform", "quadratic", not a "%s"'];
      error(fmt, convType)
  end
end

% Task 5
function fourierApprox(f,a,b,n,meth,filename)
  L = (b-a)/2;
  x  = linspace(a,b,1000);
  xx = linspace(-L,L,1000);
  partSum = zeros(1,length(x));
  
  switch meth
    case "trigonometric"
      basisK = @(k) getTrig(f,xx,k);
    case "chebyshev"
      basisK = @(k) getCheb(k);
    case "gegenbauer"
      basisK = @(k) getGegenbauer(k);
  end

  for i = 1:n
    basis = basisK(i);
    fi = trapz(f(xx).*basis(xx));
    fi = fi ./ (norm(basis(xx))^2);
    partSum = partSum + fi.*basis(x);
    
    %Solve problem with quality while parallel processing
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(x, f(x), x, partSum, 'LineWidth', 2)
    xlabel('x')
    xlim([a b])
    ylim([-2 2]) % Hardcode values to create more 'stable' gif
    drawnow
    
    gif(i, filename)
  end
end
function func = getTrig(f,xx,n)
  l = xx(end);
  an = trapz(f(xx).*cos(pi.*n.*xx./l));
  bn = trapz(f(xx).*sin(pi.*n.*xx./l));
  phi = atan2(bn,an);
  func = @(x) cos(pi*n*x/l - phi);
end
function func = getCheb(n)
  func = @(x) chebyshevT(n,x)./(1-x.^2).^(1/4);
end
function func = getGegenbauer(n)
  a = 20; 
  func = @(x) gegenbauerC(n,a,x).*(1-x.^2).^(a/2-1/4);
end

% Task 8
function drawSet(rho,N,a,b)
  phi = linspace(0, 2*pi, N);
  P = zeros(2,N);
  val = zeros(1,N);
  
  for i = 1:N
    [val(i), P(:,i)] = rho([cos(phi(i)); sin(phi(i))]);
  end
  %P = createPoly(P);
  A = zeros(2,2);
  B = zeros(2,1);
  p = zeros(2,N+1);
  for i = 1:N-1
    A(1,:) = P(:,i);
    A(2,:) = P(:,i+1);
    B(1,:) = val(i);
    B(2,:) = val(i+1);
    p(:,i) = linsolve(A,B);
  end
  A(1, :) = P(:, end);
  A(2, :) = P(:, 1);
  B(1, :) = val(end);
  B(2, :) = val(1);
  p(:, N) = linsolve(A, B);
  plot(P(1,:), P(2,:), '.');
  hold on 
  plot(p(2,:), -p(1, :), 'r-')
  xlim([a b])
  ylim([a b])
  xlabel('x')
  ylabel('y')
end
function [val,point] = supportCircle(x,a,r)
  val = dot(x,a)+r*norm(x);
  point = a + (r/norm(x)).*x;
end
function [val,point] = supportEllipse(x,p,P)
  val = dot(x,p) + sqrt(dot(x,P*x));
  point = p + (P*x)./sqrt(dot(x,P*x));
end
function [val,point] = supportSquare(x,a,r)
  val = dot(x,a) + r*sum(abs(x));
  point = a + r.*sign(x);
end
function [val,point] = supportRhombus(x,a,r)
    val = dot(x,a)+r*max(abs(x));
    delta = diff(abs(x));
    stepVec = [heaviside(-delta); heaviside(delta)];
    point = a + r.*diag(sign(x))*stepVec;
end

% Task 9
function [val,point] = supportLebesgue(f,l,x0)
  A = [];
  b = [];
  Aeq = [];
  beq = [];
  lb = [];
  ub = [];
  options = optimoptions('fmincon','Display','off');
  [x,fval] = fmincon(@(x) -dot(l,x),x0,A,b,Aeq,beq,lb,ub,...
              @(x) constraints(f,x), options);
  val = -fval;
  point = x;
end
function [c,ceq] = constraints(f,x)
  c = f(x);
  ceq = f(x);
end

% Task 10
function drawPolar(rho,N)
  a = -4;
  b = 4;
  H = [0; 0];
  len = 1;
  
  while len < N
    X = a + (b-a).*rand();
    Y = a + (b-a).*rand();
    [val,~] = rho([X; Y]);
    if abs(val - 1) <= 0.01
      H = cat(2,H,[X; Y]);
      len = len + 1;
    end
  end
  shp = alphaShape(H');
  shp.Alpha = 10;
  
  plot(shp);
  hold on
  drawSet(rho,N,a,b);
  hold off
end

% Service
function gif(idx, filename)
  [ax,rect] = getDrawArea();
  frame = getframe(ax,rect);
  im = frame2im(frame);
  [A,map] = rgb2ind(im,256);
  if idx == 1 
      imwrite(A,map,filename,'gif','Loopcount',inf,'DelayTime',0.1); 
  else 
      imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1); 
  end
end
function [ax, rect] = getDrawArea()
  ax = gca;
  ax.Units = 'pixels';
  pos = ax.Position;
  ti = ax.TightInset;
  rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
end
function demo(rho,N)
  L = length(rho);
  tiledlayout(L/2,L/2)
  for idx = 1:L
    nexttile
    drawSet(rho{idx},N,-4,4)
  end
end