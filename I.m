%% Task 1
a = -5;
b = 5;
N = 500;
%
f = @(x) sin(log(1+abs(x)) - x.^2);
nGrid = linspace(a,b,N);
[fMax, indMax] = max(f(nGrid));
[fMin, indMin] = min(f(nGrid));
fprintf('Max: %f, Min: %f\n', fMax, fMin);
%
plot(nGrid, f(nGrid))
xlabel('x')
ylabel('f(x)')
hold on
plot(nGrid(indMax), fMax, 'r*')
plot(nGrid(indMin), fMin, 'r*')
%% Task 2
input_nat = str2double(input('Enter natural: ','s'));
if ~isNatural(input_nat)
  return
end
%
vec = (1:2:input_nat);
vec = vec(~mod(vec,9));
%
mat = transpose(1:input_nat)*ones(1,input_nat);
%
LEN_C = input_nat*(input_nat+1);
shapeB = [input_nat,input_nat+1];
B = (1:LEN_C);
B = reshape(B,shapeB);
c = reshape(B,[1,LEN_C]);
D = B(:,shapeB);
%
disp('1: '); display(vec);
disp('2: '); display(mat);
disp('3: '); display(B); display(c); display(D);
%% Task 3
N = 315;
p = zeros(1,N);
A = ceil(rand(7)*N);
% Check distribution
% M = 100000000;
% for i = 1:M
%   r = ceil(rand()*N); % Uniform distribution
%   p(r) = p(r) + 1;
% end
% p = p/M;
% disp(norm(p - 1/N))
%
maxDiag = max(diag(A));
prodToSumRatio = prod(A,2) ./ sum(A,2);
maxProdToSumRatio = max(prodToSumRatio);
minProdToSumRatio = min(prodToSumRatio);
A = sortrows(A);
fprintf('Max diagonal: %d Max ratio: %d Min ratio: %d\n',...
        maxDiag, maxProdToSumRatio, minProdToSumRatio)
disp('Sorted matrix:')
display(A)
%% Task 4
N = 2;
M = 3;
X = rand(1,N);
Y = rand(1,M);
multMat = @(x,y) transpose(x).*y;
display(multMat(X,Y))
%% Task 5
input_prime = str2double(input('Enter prime: ','s'));
if ~isNatural(input_prime)
  return
elseif ~isprime(input_prime)
  disp('Error: Not a prime');
  return
end
%
A = randn(input_prime);
b = randn([input_prime, 1]);  
if abs(det(A)) <= eps
  disp('Error: Matrix is singular');
  return
end
% Check if unique solution exists(Kronecker–Capelli theorem)
if rank([A,b]) ~= rank(A)
  disp('Error: System is contradictory');
end
% Matrix method(Direct)
matrixSolution = A\b;
disp('Matrix method error:')
disp(norm(A*matrixSolution - b));
% Least squares method(Iterative)
lsqrSolution = lsqr(A,b,1e-14,1000);
disp('Least squares method error:')
disp(norm(A*lsqrSolution - b));
clear;
%% Task 6
N = 2;
M = 3;
vecA = rand(1,N);
vecB = rand(M,1);
max(max(vecA-vecB,vecB-vecA),[],'all')
%% Task 7
N = 3;
K = 2;
matA = rand(N,K);
matProd = matA * matA';
matDist = sqrt(ones(N,1)*diag(matProd)'+diag(matProd)*ones(1,N)-2.*matProd);
display(matDist);
%% Task 8
input_nat = str2double(input('Enter natural: ','s'));
if ~isNatural(input_nat)
  return
end
matA = dec2bin(0:1:2^input_nat-1) - '0';
display(matA);
%% Task 9
N = 100;
timesMat = zeros(2,N);
for i = 1:N
  A = rand(i);
  B = rand(i);
  timesMat(1,i) = compTimeMedian(@(A,B) myMultiply(A,B),A,B);
  timesMat(2,i) = compTimeMedian(@(A,B) mtimes(A,B),A,B);
end
plot((1:N),timesMat(1,:),(1:N),timesMat(2,:),'LineWidth',3)
legend('My multiplication', 'Default multiplication')
%% Task 10
X = [NaN 1 2; NaN 0 6; 1 5 NaN];
display(mean(X,'omitnan'))
%% Task 11
N = 1000;
a = ceil(N*rand);
sigma = ceil(N*rand./4);
vec = sigma.*randn(1,N) + a;
count3sigma = sum(vec >= a-3.*sigma & vec <= a+3.*sigma);
disp(count3sigma./N);
%% Task 12
a = -4*pi;
b = 4*pi;
h = 0.001;
%
n = ceil((b-a)/h);
X = linspace(a,b,n);
%
f = @(x) sin(x)./x;
intTrapz = @(x) trapz(x, f(x));
intRect  = @(x) rectangles(x, f(x));
intSimps = @(x) simpson(x, f(x));
Xi = @(H, i) (a:H:X(i));
%
Y = zeros(3,n);
timeVec = zeros(3,n);
for i = 2:n
  tic
  Y(1,i) = intTrapz(Xi(h,i));
  timeVec(1,i) = toc;
  tic
  Y(2,i) = intRect(Xi(h,i));
  timeVec(2,i) = toc;
  tic
  Y(3,i) = intSimps(Xi(h,i));
  timeVec(3,i) = toc;
end
%
timeVec = 1000000*timeVec;
fprintf('Time: Trapz=%.2fμs Rectangular=%.2fμs Simpson=%.2fμs\n',...
        median(timeVec(1,:)),median(timeVec(2,:)),median(timeVec(3,:)));
%
convF = @(F,H,i) F(Xi(H,i)) - F(Xi(H/2,i));
Hi = linspace(eps, 0.1, n);
convRate = zeros(3,n);
for j = 2:n
  convRate(1,j) = convF(intTrapz,Hi(j),j);
  convRate(2,j) = convF(intRect, Hi(j),j);
  convRate(3,j) = convF(intSimps,Hi(j),j);
end
%
tiledlayout(2,1)
% Top plot
nexttile
plot(X,Y(1,:),X,Y(2,:),'--',X,Y(3,:),':','LineWidth',3);
xlabel('x')
ylabel('f(x)')
legend('Trapz','Rectangular','Simpson')
% Bottom plot
nexttile
plot(Hi, convRate(1,:), ...
     Hi, convRate(2,:), '--',...
     Hi, convRate(3,:), ':', 'LineWidth',3)
xlabel('h')
ylabel('Convolution Rate(h)')
legend('Trapz','Rectangular','Simpson')
%% Task 13
H = logspace(-8,-1,1000);
X = 1;
%
f = @(x) erf(x);
fDiff = @(x) 2/sqrt(pi)*exp(-x^2);
fRightDiff = @(x, h) (f(x+h) - f(x))./h;
fCentDiff  = @(x, h) (f(x+h) - f(x-h))./(2.*h);
absDiff = @(F) abs(fDiff(X) - F(X,H));
loglog(H, absDiff(fRightDiff), ... 
       H, absDiff(fCentDiff), '--', 'LineWidth',3)
xlabel('h')
ylabel('Abs Diff')
legend('Right Diff', 'Centr Diff')
%% Functions
function res = isNatural(n)
  res = 0;
  if isnan(n)
    disp('Error: Not a number');
    return
  elseif n <= 0 
    disp('Error: Non-positive');
    return
  elseif mod(n, 1) ~= 0
    disp('Error: Non an integer');
    return
  end
  res = 1;
end
%
function C = myMultiply(A,B)
  szA = size(A);
  szB = size(B);
  C = zeros(szA(1),szB(2));
  for i = 1:szA(1)
    for j = 1:szB(2)
      C(i,j) = A(i,:)*B(:,j);
    end
  end
end
%
function t = compTimeMedian(F,A,B)
  numTries = 10;
  res = zeros(numTries, 1);
  for i = 1:numTries
    tic();
    F(A,B);
    res(i) = toc();
  end
  t = median(res);
end
%
function res = rectangles(X,Y)
  n = length(X);
  h = (X(n)-X(1))/n;
  res = h*(Y(1)/2 + sum(Y(2:end-1)) + Y(end)/2);
end
%
function res = simpson(X,Y)
  n = length(X);
  h = (X(n)-X(1))/n;
  res = Y(1)+2*sum(Y(3:2:end-2))+4*sum(Y(2:2:end))+Y(end);
  res = h/3*res;
end