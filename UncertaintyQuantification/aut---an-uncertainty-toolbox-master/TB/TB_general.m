% test bench --> should at least not raise error

a=unc(1,1);
b=unc(1,1);
s=a+b;
m=a*b;
d=a/a;
r=a-b;

v=[a b a b];

v1=[a s s a];

%A3=unc([1 2 3;1 2 3;3 2 1])+a;
A3=unc([1 2 3;1 2 3;3 2 1],[1 2 3;1 2 3;3 2 1]*0.1)+a;

A=sqrt(A3);

A3(1,1)=A3(1,1)+1;
A3=A3+a;
A3=A3-a;

A3*v(1:3)';

A3 = [unc(1,0.1),unc(1,0.1),unc(1,0.1);...
      unc(2,0.1),unc(3,0.1),unc(4,0.1);...
      unc(3,0.1),unc(4,0.1),unc(5,0.1)];


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

A3.^4;

A3./A3;

3*A;
3+A;
2+v;


A3(1,1)=A3(1,1)+1;
inv(A3);
inv(A3');
A3\v(1:3)';

sqrt(a);
sqrt(v);
sqrt(A);







