%% Task 1
N = 100;
X = linspace(0, 2, N); 
f = @(x) x.^3 + x - 1;
g = @(x) sin(x);
h = @(x) f(x) - g(x);
plot(X, f(X), X, g(X), 'LineWidth', 2)
[Pinput , ~] = ginput;
for i = 1:size(Pinput)
    [xi, hval] = fzero(h, Pinput(i));
    fprintf("X: %e, \x0394: %e\n", xi, abs(hval))
end
fprintf('\n') % For multiple executions
%% Task 2
N = 100;
alpha = 0.001;
X = linspace(-alpha, alpha, N);
f = @(x) 0.*(abs(x) <= eps) + x.*cos(log(abs(x))).*(abs(x) >= eps);
plot(X, f(X), 'LineWidth', 1.5)
Y = arrayfun(@(x) fzero(f,x), X);
hold on
plot(X, Y, 'LineWidth', 1.5)
hold off
legend('xcos(ln(|x|))', 'Root approximation')
%% Task 3
N = 1000;
A = [pi exp(1); exp(1) pi];
I = eye(2);

TaylorApprox = I;
for i = 1:N
    I = I*A/i;
    TaylorApprox = TaylorApprox + I;
end

B = [A, zeros(2); zeros(2), A];
[~, OdeApprox] = ode45(@(t,y) B * y, linspace(0, 1, N), [1, 0, 0, 1]);
OdeApprox = reshape(OdeApprox(end,:), 2, 2);

MatlabApprox = expm(A);

absErrorTaylor = abs(MatlabApprox - TaylorApprox);
relErrorTaylor = max(absErrorTaylor(:) ./ abs(MatlabApprox(:)));
absErrorOde = abs(MatlabApprox - OdeApprox);
relErrorOde = max(absErrorOde(:) ./ abs(MatlabApprox(:)));

fprintf("Relative errors(compared to Matlab expm()) with N=%d\n", N)
fprintf("Taylor: %e, Ode(ode45): %e\n\n",relErrorTaylor,relErrorOde)
%% Task 4
cla reset;

N = 1000;
alpha = 0.1;
rectWidth = 2;
rectHeight = 3;
x0 = 1;
y0 = 1;
x1 = 3;
y1 = 1;

dxdt = @(t, x) [x(3); x(4); -alpha*x(1); -alpha*x(2)];    % Standart subsitution
[~, X] = ode45(dxdt, linspace(0, 10, N), [x0; y0; x1; y1]);
X = computeCollisions(X, rectWidth/2, rectHeight/2);      % Using symmetry(/2)

plotRectangle(rectWidth/2, rectHeight/2)
hold on
comet(X(:,1), X(:,2))
%% Task 5
N = 128;
T = linspace(0, 10, N);
[~, X] = ode45(@gravityDiff, T, [0; 0; 10; 0; 0; 10; 0; -10]);

minX = min([X(:,1); X(:,3)]);
maxX = max([X(:,1); X(:,3)]);
minY = min([X(:,2); X(:,4)]);
maxY = max([X(:,2); X(:,4)]);
mov(1:N) = struct('cdata', [],'colormap', []);
for i = 1:N
    plot(X(1:i,1), X(1:i,2), X(1:i,3), X(1:i,4))
    axis([minX maxX minY maxY])
    mov(i) = getframe();
end
%% Task 6
SavePrivateRyan(1);
%% Task 7 (I)
N = 100;
a = -1;
b = 1;
V = @(x, y) x.^2 + y.^2;                        % Lyapunov function
Phi = linspace(0, 2*pi, 10);                    % N = 10 guarantee optimal computation time

x0 = [cos(Phi); sin(Phi)] / 2;
x = linspace(a, b, N);
y = linspace(a, b, N);
T = linspace(0, 2, N);

[X, Y] = meshgrid(x, y);
contour(X, Y, V(X, Y), 20, '--');
hold on;
dydt = @(t, y) [y(1)^3 - y(2); y(1) + y(2)^3];
for i = 1:size(x0, 2)
    [T, Y] = ode45(dydt, T, x0(1:2, i));
    maxV = max(V(Y(:,1), Y(:,2)));
    for j = 1:size(Y, 1)-2
        quiver(Y(j,1), Y(j,2), Y(j+2,1)-Y(j,1), Y(j+2,2)-Y(j,2), ...
                'Color', [V(Y(j,1), Y(j,2))/maxV, 0, 0]);
    end
end
xlim([a b])
ylim([a b])
%% Task 7 (II)
N = 10;
a = -1;
b = 1;
V = @(x, y) x.^2 + y.^4;          % Lyapunov function

Phi = linspace(0, 2*pi, N);
x0 = [cos(Phi); sin(Phi)] ;
x = linspace(a, b, 100);
y = linspace(a, b, 100);
T = linspace(0, 2, N);

[X, Y] = meshgrid(x, y);
contour(X, Y, V(X, Y), 20, '--');
hold on;
dydt = @(t, y) [y(2)^3 - y(1)^5; -y(1) - y(2)^3 + y(2)^5];
for i = 1:size(x0, 2)
    [~, Y] = ode45(dydt, T, x0(1:2, i));
    maxV = max(V(Y(:,1), Y(:,2)));
    for j = 1:size(Y, 1)-2
        quiver(Y(j,1), Y(j,2), Y(j+2,1)-Y(j,1), Y(j+2,2)-Y(j,2), ...
                'Color', [V(Y(j,1), Y(j,2))/maxV, 0, 0]);
    end
end
ylim([a b])
xlim([a b])
%% Task 8
N = 100;
syms y(t)

cond1 = y(0) == 0;
cond2 = y(1) == -1;
conds = [cond1 cond2];
ode = diff(y,t,2) - y(t) == 2*t;
ySol(t) = dsolve(ode, conds);
analyticalSol = matlabFunction(ySol);

dydx = @(x,y) [y(2), y(1) + 2*x];
res = @(ya,yb) [ya(1), yb(1)+1];
guess = @(x) [sin(x), cos(x)];
xmesh = linspace(0,1,N);
solinit = bvpinit(xmesh, guess);
sol = bvp4c(dydx, res, solinit);

errorL2 = sqrt( trapz((sol.y(1,:) - analyticalSol(sol.x)).^2) );
errorC  = max( abs(sol.y(1,:) - analyticalSol(sol.x)) );
fprintf("L2: %e, C: %e\n\n", errorL2, errorC)
%% Task 9
f = @(x) sin(x(1)) / x(1);
g = @(x) x(1).^2 + x(1)*x(2) + x(2).^2 - 3*x(1);
h = @(x) (x(1) + x(2)).^2 + (x(3)-1).^2 + 1;

fmin(f, 1)
fmin(g,[2 -1])
fmin(h,[1 2 3])
disp("____________________")
%% Task 10
f1 = @(t) exp(-3.*abs(t)).*(sin(t).^3);
f2 = @(t) (sin(t) - t.*cos(t)) ./ t.^2;
f3 = @(t) cos(t)/(1+abs(t.^3));
f4 = @(t) exp(-5.*t.^8).*sin(t+t.^3);
f2FT = @(t) -1i*(pi/2).*t.*(sign(1-t)+sign(1+t));

X  = linspace(-4*pi,4*pi,1000);
F = f1FT(X);
F3 = repmat(F, 1, 3);
X3 = linspace(-35,35,3000);

fighandle1 = figure('Name','f1');
res = plotFT(fighandle1, f1, @f1FT, 0.001,[-100 100],[-35 35]);
% fighandle1 = figure('Name','f2');
% res = plotFT(fighandle1, f2, f2FT, 0.05, [-50 50], [-2 2]);

% fighandle1 = figure('Name','f3');
% % res = plotFT(fighandle1, f3, [], 0.1, [-1 1], [-5 5]);

% fighandle1 = figure('Name','f4');
% res = plotFT(fighandle1, f4, [], 0.001, [-100 100], [-20 20]);
saveas(gcf,'Aliasing-2.png')
%%
fighandle = figure('Name','func2');
subplot(2, 1, 1);
ax = gca;
ax.XLim = [-4 3];
subplot(2, 1, 2);
ax = gca;
ax.XLim = [-4 3];
SPlotInfo = struct('realXLim', 1, 'imagXLim', 1);
set(fighandle,'UserData',SPlotInfo);
res = plotFT(fighandle, func2, ftfunc2, 0.01, [-10 20], [-5 5]);
%%
fighandle = figure('Name','func2');
SPlotInfo = struct('outLimVec', [0 5], 'inpLimVec', [-10 10]);
set(fighandle,'UserData',SPlotInfo);
res = plotFT(fighandle, func2, ftfunc2, 0.01, [-10 20], [-3 5]);
%% Functions
function plotRectangle(rectWidth, rectHeight)
  ident = 1;
  argLim = max(rectWidth + ident,rectHeight + ident);

  hold on
  plot([-rectWidth rectWidth], [rectHeight rectHeight], 'black')
  plot([-rectWidth rectWidth], [-rectHeight -rectHeight], 'black') 
  plot([rectWidth rectWidth],  [-rectHeight rectHeight], 'black')
  plot([-rectWidth -rectWidth],[-rectHeight rectHeight], 'black')
  axis([-argLim argLim -argLim argLim])
  xlabel('x');
  ylabel('y');
  hold off
end
function Y = computeCollisions(X, rectWidth, rectHeight)
  N = size(X,1);
  Y = X;
  collisionX = 1;
  collisionY = 1;
  for i = 1:N
      if abs(Y(i, 1)) > rectWidth && collisionX
          Y(i:end, 1) = 2*Y(i,1)-Y(i:end, 1);
          collisionX = 0;
      else
          collisionX = 1;
      end
      if abs(Y(i, 2)) > rectHeight && collisionY
          Y(i:end, 2) = 2*Y(i,2)-Y(i:end, 2);
          collisionY = 0;
      else
          collisionY = 1;
      end
  end
end

function res = gravityDiff(~,x)
  m1 = 1;
  m2 = 10;
%     m1 = 1.989 * 10^30;
%     m2 = 5.972 * 10^24;
%     G = 6.67384*10^(-11);
  G = 300;

  res = zeros(8, 1);
  res(1:2) = x(5:6);
  res(3:4) = x(7:8);
  res(5) = G*m2/(norm(x(1:2) - x(3:4))^3) * (x(3) - x(1));
  res(6) = G*m2/(norm(x(1:2) - x(3:4))^3) * (x(4) - x(2));
  res(7) = G*m1/(norm(x(1:2) - x(3:4))^3) * (x(1) - x(3));
  res(8) = G*m1/(norm(x(1:2) - x(3:4))^3) * (x(2) - x(4));
end

function fmin(f, x0)
  N = 1000;
  a = 5;
  argNum = size(x0, 2);
  if argNum == 2
      x = linspace(-a, a, N);
      y = linspace(-a, a, N);
      [X,Y] = meshgrid(x, y);
      Z = zeros(N, N);
      for i = 1:N
          for j = 1:N
              Z(i,j) = f([x(i) y(j)]);
          end
      end
      contour(X, Y, Z, 10);
      hold on
      prevX0 = x0;
  end
  opts = optimoptions('fminunc', 'Display', 'off');
  for i = 1:argNum
      g = @(x) f([x0(1:(i-1)) x x0(i+1:argNum)]);
      x0(i) = fminunc(g, x0(i), opts);
      if argNum == 2
          plot([prevX0(1) x0(1)], [prevX0(2) x0(2)], '--')
      end
      prevX0 = x0;
  end
  fprintf("minimum: %f\n", f(x0))
  if argNum == 1
    fprintf("fminbnd: %f\n", f(fminbnd(f, -a, a)))
  end
end

function GenerateTable(n,m)
  P = rand(n,m);
  save('table.mat', 'P')
  tiledlayout(2,1)
  nexttile
  image(P,'CDataMapping','scaled')
  colorbar
  disp(P)
end                                                             
function mat = neigh(mat, i, j, var)
  [n,m] = size(mat);
  if var
      if i + 1 <= n && j + 1 <= m
          mat(i+1, j+1) = 1;
      end
      if i + 1 <= n && j - 1 > 0
          mat(i+1, j-1) = 1;
      end
      if i - 1 > 0 && j + 1 <= m
          mat(i-1, j+1) = 1;
      end
      if i - 1 > 0 && j - 1 > 0
          mat(i-1, j-1) = 1;
      end
  end
  if i + 1 <= n
      mat(i+1, j) = 1;
  end
  if i - 1 > 0
      mat(i-1, j) = 1;
  end
  if j + 1 <= m
      mat(i, j+1) = 1;
  end
  if j - 1 > 0
      mat(i, j-1) = 1;
  end
end     
function val = min_neigh(fld, P, i, j, var)
  [n, m] = size(P);
  val = -Inf;
  if var
      if i + 1 <= n && j + 1 <= m && P(i, j) * fld(i+1, j+1) > val
          val = P(i, j) * fld(i+1, j+1);
      end
      if i + 1 <= n && j - 1 > 0 && P(i, j) * fld(i+1, j-1) > val
          val = P(i, j) * fld(i+1, j-1);
      end
      if i - 1 > 0 && j + 1 <= m && P(i, j) * fld(i-1, j+1) > val
          val = P(i, j) * fld(i-1, j+1);
      end
      if i - 1 > 0 && j - 1 > 0 && P(i, j) * fld(i-1, j-1) > val
          val = P(i, j) * fld(i-1, j-1);
      end
  end
  if i + 1 <= n && P(i, j) * fld(i+1, j) > val
      val = P(i, j) * fld(i+1, j);
  end
  if i - 1 > 0 && P(i, j) * fld(i-1, j) > val
      val = P(i, j) * fld(i-1, j);
  end
  if j + 1 <= m && P(i, j) * fld(i, j+1) > val
      val = P(i, j) * fld(i, j+1);
  end
  if j - 1 > 0 && P(i, j) * fld(i, j-1) > val
      val = P(i, j) * fld(i, j-1);
  end
end
function dir = maxval(fld, var, cur)
  maxv = -Inf;
  i = cur(1);
  j = cur(2);
  [n,m] = size(fld);
  if var
      if i + 1 <= n && j + 1 <= m && maxv < fld(i+1, j+1)
          maxv = fld(i+1, j+1);
          dir = [i+1, j+1];
      end
      if i + 1 <= n && j - 1 > 0 && maxv < fld(i+1, j-1)
          maxv = fld(i+1, j-1);
          dir = [i+1, j-1];
      end
      if i - 1 > 0 && j + 1 <= m && maxv < fld(i-1, j+1)
          maxv = fld(i-1, j+1);
          dir = [i-1, j+1];
      end
      if i - 1 > 0 && j - 1 > 0 && maxv < fld(i-1, j-1)
          maxv = fld(i-1, j-1);
          dir = [i-1, j-1];
      end
  end
  if i + 1 <= n && maxv < fld(i+1, j)
      maxv = fld(i+1, j);
      dir = [i+1, j];
  end
  if i - 1 > 0 && maxv < fld(i-1, j)
      maxv = fld(i-1, j);
      dir = [i-1, j];
  end
  if j + 1 <= m && maxv < fld(i, j+1)
      maxv = fld(i, j+1);
      dir = [i, j+1];
  end
  if j - 1 > 0 && maxv < fld(i, j-1)
      dir = [i, j-1];
  end
end
function SavePrivateRyan(var)
  n = 5;
  m = 7;
  GenerateTable(n,m);
  load('table.mat')
  fld = -inf*ones(n, m);
  fld(n,m) = 1;
  change = zeros(n,m);
  change(n-1, m) = 1;
  change(n, m-1) = 1;
  while 1
      new_change = zeros(n,m);
      for i = n:-1:1
          for j = m:-1:1
              if change(i,j)
                  val = min_neigh(fld, P, i, j, var);
                  if val > fld(i,j)
                      new_change = neigh(new_change, i, j, var);
                      fld(i,j) = val;
                  end
              end
          end
      end
      if min(new_change == change)
          break
      end
      change = new_change; 
  end
  
  disp(fld)
  nexttile
  image(fld,'CDataMapping','scaled')
  colorbar
  disp(fld)
  hold on
   for i = 1:n
       for j = 0:(m-1)
           text(0.6 + j,i - 0.2, string(fld(i,j + 1)))
       end
   end
  direct = [1, 1];
  while min(direct == [n, m]) ~= 1
      new_direct = maxval(fld, var, direct);
      plot([direct(2), new_direct(2)], [direct(1), new_direct(1)], 'r');
      direct = new_direct;
  end
end

function res = plotFT(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
  figure(hFigure);
  if isfield(hFigure.UserData, 'inpLimVec')
      a = hFigure.UserData.inpLimVec(1);
      b = hFigure.UserData.inpLimVec(2);
  else
      a = inpLimVec(1);
      b = inpLimVec(2);
  end
  if isfield(hFigure.UserData, 'outLimVec')
      c = hFigure.UserData.outLimVec(1);
      d = hFigure.UserData.outLimVec(2);
  elseif length(outLimVec) > 1 
      c = outLimVec(1);
      d = outLimVec(2);
  else
      c = -1/(2*step);
      d = 1/(2*step);
  end
  info = get(hFigure, 'UserData');
  if ~isempty(info) && isfield(hFigure.UserData, 'ax_handles')
      delete(info.ax_handles);
  else
      info = struct('first_axises', [], 'second_axises', []);
  end
  T = b - a;
  N = round(T / step);
  step = T / N;
  x = a:step:b;
  y = fHandle(x);
  y(isnan(y)) = 0;
  n = 0;
  if a > 0
      n = round(a/T);
  elseif b < 0
      n = -round(b/T);
  end
  border = n*T;
  ind = find(x <= border, 1, 'last');
  Y(1:N+1-ind) = y(ind+1:N+1);
  Y(N+1-ind+1:N+1) = y(1:ind);
  Yft = fft(Y)*T/(N+1);
  stepft = 2*pi/T;
  Tft = stepft*(N+1);
  lbord = floor(c / Tft)*Tft;
  rbord = ceil(d / Tft)*Tft;
  counter = ceil(d / Tft) - floor(c / Tft);

  Xft = lbord:stepft:rbord - stepft;
  Yft = repmat(Yft(1:end), 1, counter);
  subplot(2, 1, 1);
  hold on;
  plot(Xft(1:end), real(Yft), 'b','LineWidth',1.5);
  legend('Fourier transform approximation');
  if ~isempty(fFTHandle)
      res_f = fFTHandle(linspace(c, d, 1000));
      plot(linspace(c, d, 1000), real(res_f), 'r', 'DisplayName', 'Analytical Fourier transform','LineWidth',1.5);
      if isfield(hFigure.UserData, 'realXLim') 
          ylim([min(real(Yft))*1.2 max(real(Yft))*1.2])
      else
          axis([c d min([real(Yft) real(res_f)])   max([real(Yft) real(res_f)])]);
      end
  elseif isfield(hFigure.UserData, 'realXLim') 
      ylim([min(real(Yft))*1.2 max(real(Yft))*1.2])
  else
      axis([c d min(real(Yft)) max(real(Yft))]);
  end
  ylabel('Re F(\lambda)');
  xlabel('\lambda');
  title('Real axis');
  info.first_axises = gca;  
  subplot(2, 1, 2);
  hold on;
  plot(Xft(1:end), imag(Yft), 'b','LineWidth',1.5);
  legend('Fourier transform approximation');
  if ~isempty(fFTHandle)
      res_f = fFTHandle(linspace(c, d, 1000));
      plot(linspace(c, d, 1000), imag(res_f), 'r', 'DisplayName', 'Analytical Fourier transform','LineWidth',1.5);
      if isfield(hFigure.UserData, 'imagXLim') 
          ylim([min(imag(Yft))*1.2 max(imag(Yft))*1.2])
      else
          axis([c d min([imag(Yft) imag(res_f)])   max([imag(Yft) imag(res_f)])]);
      end
  elseif isfield(hFigure.UserData, 'imagXLim') 
      ylim([min(imag(Yft))*1.2 max(imag(Yft))*1.2])
  else
      axis([c d min(imag(Yft)) max(imag(Yft))]);
  end
  xlim([c d])
  title('Imaginary axis');
  ylabel('Im F(\lambda)');
  xlabel('\lambda');
  info.second_axises = gca;
  set(hFigure,'UserData',info);
  res.inpLimVec = [a b];
  res.outLimVec = [c d];
end
function res = f1FT(w)
  a = 3.75e-1i./(w.*1i+-3.0-1i);
  b = 3.75e-1i./(w.*1i+-3.0+1i);
  c = 3.75e-1i./(w.*1i+3.0-1i);
  d = 3.75e-1i./(w.*1i+3.0+1i);
  e = 1.25e-1i./(w.*1i+-3.0-3.0i);
  f = 1.25e-1i./(w.*1i+-3.0+3.0i);
  g = 1.25e-1i./(w.*1i+3.0-3.0i);
  h = 1.25e-1i./(w.*1i+3.0+3.0i);
  res = a-b-c+d-e+f+g-h;
end