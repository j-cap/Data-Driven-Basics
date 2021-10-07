% defining uncertain objects 
a=unc_t(1,1);
b=unc_t(1,1);
% basic operation with real numbers 
s=a+b; 
m=a*b;
d=a/a;
r=a-b;
% defining uncertain vectors and matrices 
v=[a b a b];

v1=[a s s a];

A3=unc_t([1 2 3;1 2 3;3 2 1],[1 2 3;1 2 3;3 2 1]*0.1)+a;

% operation with matrices 

A=sqrt(A3);


A3(1,1)=A3(1,1)+1;
A3=A3+a;
A3=A3-a;

A3+v(1:3);

A3*v(1:3)';

A3 = [unc_t(1,0.01),unc_t(1,0.01),unc_t(1,0.01);...
      unc_t(2,0.01),unc_t(3,0.01),unc_t(4,0.01);...
      unc_t(3,0.01),unc_t(4,0.01),unc_t(6,0.01)];


Res = A3-A3;
if (any(any(gmv(Res)~=zeros(3,3))) || any(any(gmu(Res)~=zeros(3,3))))
    ErrorFlag = 1;
    disp ('Subtraction Error! (A3-A3)');
end

Res = A3./A3;
if size(Res,1)~=3 || size(Res,2)~=3
    ErrorFlag = 1;
    disp ('Division Error! (A3./A3) Matrix dim mismatch');    
else
    if (any(any(gmv(Res)~=ones(3,3))) || any(any(gmu(Res)~=zeros(3,3))))
        ErrorFlag = 1;
        disp ('Division Error! (A3./A3)');
    end
end

A3*A3;

sum([A3,A3]);
cumsum([A3,A3]);

sum([A3,A3]');
cumsum([A3,A3]');

sin(A3)./A3;

(A3+1)*3./a;

tic 
A3.^4
A3.^v(1:3)
toc

A3./A3;

3*A;
3+A;
2+v;

%% test inverse function 
A = rand(3,3)
Au = rand(3,3);
A_mc = unc_t(A,Au);  % MC approach
A_lin = unc_t(A,Au); % linearization approach 

A_I = inv (A);
A_mcI = inv(A_mc);
A_linI = inv(A_lin);

if (any(any((gmv(A_mcI)- A_I)~=zeros(3,3))))
    ErrorFlag = 1;
    disp ('Check result of inverse matrix for MC approach!');
elseif (any(any((gmv(A_linI)- A_I)~=zeros(3,3))))
    ErrorFlag = 1;
    disp ('Check result of inverse matrix for lineatization approach approach!');
end

A3\v(1:3)';

sqrt(a);
sqrt(v);
sqrt(A);
