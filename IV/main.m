%% Task 1
N = 3;
Av = complex(rand(1, N), rand(1,N));
Bv = complex(rand(1, N), rand(1,N));
Cv = complex(rand(1, N), rand(1,N));

Pv = [Av; Bv; Cv];
Xv = zeros(N, 2);

for i = 1:N
  Xv(i,:) = roots(Pv(:,i))';
end
[X1, X2,D] = quadsolve(complex(Av),complex(Bv),complex(Cv));
disp(Xv')
disp([X2; X1])
disp('_________________________________________________________')
%% Task 2
N = 3;
A = rand(N);
I = inv(A);
Imat = inv_matlab(A);
Ic = inv_c(A);

absErrorC = abs(Ic-I);
relErrorC = max(absErrorC(:) ./ abs(Ic(:)));
absErrorMatlab = abs(Imat-I);
relErrorMatlab = max(absErrorMatlab(:) ./ abs(Imat(:)));
fprintf('C error: %d, Matlab error: %d\n', relErrorC, relErrorMatlab);
%% Task 3
N = 50;
Ngrid = 1:N;
ErrLin = zeros(1,N);
ErrMat = zeros(1,N);
ErrC   = zeros(1,N);

for i = Ngrid
  A = rand(i);
  Itrue = inv(A);
  Itest = {linsolve(A,eye(i)), inv_matlab(A), inv_c(A)};
  absError = cellfun(@(X) abs(X-Itrue), Itest, 'UniformOutput', false);
  relError = cellfun(@(X,Err) max(Err(:) ./ abs(X(:))), Itest, absError);
  
  ErrLin(1,i) = relError(1);
  ErrMat(1,i) = relError(2);
  ErrC(1,i)   = relError(3);
end
plot(Ngrid, ErrLin, Ngrid, ErrMat, Ngrid, ErrC, 'LineWidth',1.5)
legend('Linsolve','Matlab','C')
%% Task 4
N = 50;
Ngrid = 1:N;
TimeInv = zeros(1,N);
TimeLin = zeros(1,N);
TimeMat = zeros(1,N);
TimeC   = zeros(1,N);

for i = Ngrid
  A = rand(i);

  TimeInv(1,i) = exec_time(@inv,A);
  TimeLin(1,i) = exec_time(@(X) linsolve(X,eye(i)),A);
  TimeMat(1,i) = exec_time(@inv_matlab,A);
  TimeC(1,i)   = exec_time(@inv_c,A);
end
plot(Ngrid,TimeInv,Ngrid,TimeLin,Ngrid,TimeMat,Ngrid,TimeC,'LineWidth',1.5)
legend('Inverse','Linsolve','Matlab','C')
%% Task 5
N = 20;
Ngrid = 1:N;

ss = {'inv', 'linsolve', 'inv\_matlab', 'inv\_c'};    % Add backslash for legend support
regs = cellfun(@(s) fitlm(Ngrid, Time(N,s), 'quadratic'), ...
                ss,'UniformOutput', false);
cla reset;  % Clear previous screen
hold on
for i = 1:4 
  plot(Ngrid, regs{i}.Fitted, 'DisplayName', ss{i}, 'LineWidth', 2);
end
hold off
legend show
%% Task 6
M = 10;
N = 12;
mu = 0;
fHandle   = @(x,y) x+y;
xiHandle  = @(x) x.^2 - x;
etaHandle = @(y) 2*(y.^2 - y);

solveDirichlet(fHandle, xiHandle, etaHandle, mu, M, N)
%% Test
N = 100;
X = linspace(0,1,N);
K = 1:N;
f = @(x) -2.*x.*sin(x);
mu = 0.1;
uZero = 1;

syms u(t)
cond1 = u(0) == uZero;
cond2 = u(1) == uZero;
conds = [cond1 cond2];
ode = diff(u,t,2) - mu * u(t) == f(t);
uSol(t) = dsolve(ode,conds);
uF = matlabFunction(uSol);
U = uF(X);

D = -4*sin(2*pi*K/N)*N^2 - mu;
Dinv = 1./D;
F = f(X);
F(1) = 0;
Bint = ifft(F);
F0 = (uZero - Bint*Dinv')/(sum(Dinv)/N);
F(1) = F0;
F(end) = F0;
B = ifft(F);
A = B .* Dinv;
Ufft = fft(A);

tiledlayout(2,1)
nexttile
plot(X,real(Ufft))
nexttile
plot(X,U)
%%
a = 0;
b = 1;
N = 4096;
mu = 0.1;
t = linspace(a,b,N);
f = -2.*t.*sin(t);
fftx = fft(f);
k = (2*pi/(b-a))*[0:N/2-1, 0, -N/2+1:-1];% "swapping"
dffft = -fftx./(k.^2+mu);% multiplication by ik
df2 = ifft(dffft);
plot(t,real(df2),'.')
%% Functions
function I = inv_matlab(A)
  [N, M] = size(A);
  if M ~= N
    error('Square matrices only') 
  end
  
  I = eye(N);
  M = [A, I];
  % Forward alg(below main diag)
  for k = 1:N
    for i = 1:2*N
      M(k,i) = M(k,i) / A(k,k);
    end
    for i = k+1:N
      coeff = M(i,k) / M(k,k);
      for j = 1:2*N
        M(i,j) = M(i,j) - M(k,j) * coeff;
      end
    end
    for i = 1:N
      for j = 1:N
        A(i,j) = M(i,j);
      end
    end
  end
  % Backwards alg(above main diag)
  for k = N:-1:1
    for i = 2*N:-1:1
      M(k,i) = M(k,i) / A(k,k);
    end
    for i = k-1:-1:1
      coeff = M(i,k) / M(k,k);
      for j = 2*N:-1:1
        M(i,j) = M(i,j) - M(k,j) * coeff;
      end
    end
  end

  for i = 1:N
    for j = 1:N
      I(i,j) = M(i,N+j);
    end
  end
end

function t = exec_time(f, arg)
  N = 100;
  tSum = 0;
  for i = 1:N
    tic;
    f(arg);
    tSum = tSum + toc;
  end
  t = tSum/N;
end

function t = Time(n,s)
  Ngrid = 1:n;
  TimeInv = zeros(1,n);
  TimeLin = zeros(1,n);
  TimeMat = zeros(1,n);
  TimeC   = zeros(1,n);
  
  for i = Ngrid
    A = rand(i);

    TimeInv(1,i) = exec_time(@inv,A);
    TimeLin(1,i) = exec_time(@(X) linsolve(X,eye(i)),A);
    TimeMat(1,i) = exec_time(@inv_matlab,A);
    TimeC(1,i)   = exec_time(@inv_c,A);
  end
  
  switch s
    case 'inv'
      t = TimeInv;
    case 'linsolve'
      t = TimeLin;
    case 'inv\_matlab'
      t = TimeMat;
    case 'inv\_c'
      t = TimeC;
    otherwise
      error('No such method')
  end
end

function res = solveDirichlet(fHandle,xiHandle,etaHandle,mu,M,N)
  x = linspace(0,1,M);
  y = linspace(0,1,N);
  
  [X,Y] = meshgrid(x,y);
  Fsp = fHandle(X,Y); 
  Xi = xiHandle(x);
  Eta = etaHandle(y);
  
  B = ifft2(F);
end

function res = uAnalytical(xMat,yMat,u1Zero,u2Zero,mu)
  f1 = @(x) -2*x*sin(x);
  f2 = @(y) (4+y)*exp(-2*y);
  
  syms u1(t) u2(t)
  cond1 = u1(0) == u1Zero;
  cond2 = u1(1) == u1Zero;
  cond3 = u2(0) == u2Zero;
  cond4 = u2(1) == u2Zero;
  conds1 = [cond1 cond2];
  conds2 = [cond3 cond4];
  ode1 = diff(u1,t,2) - mu * u1(t) == f1(t);
  ode2 = diff(u2,t,2) - mu * u2(t) == f2(t);
  u1Sol(t) = dsolve(ode1,conds1);
  u2Sol(t) = dsolve(ode2,conds2);
  
  u1F = matlabFunction(u1Sol);
  u2F = matlabFunction(u2Sol);
  u = @(x,y) u1F(x) + u2F(y);
  res = u(xMat,yMat);
end

function res = fGiven(x,y)
  res = -2*x*sin(x)+(4+y)*e^(-2*y);
end