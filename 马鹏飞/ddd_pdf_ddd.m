clc
clear
alpha=1.0;
beta=0.8;
eta=0.6;
sig2=0.1;

%%%%%%%%%%%%%%%%%%%%%%%
dx=0.02;dy=0.02;
lx=-5;rx=5;
ly=-5;ry=5;
x=lx:dx:rx;
y=ly:dy:ry;


test1=0;
for i=1:length(x)
    m=0;
    for j=1:length(y)
        % p=exp(mu/(2*D)*((sqrt(x(i)^2+y(j)^2))^2-1/8*(sqrt(x(i)^2+y(j)^2))^4+1/24*alpha*(sqrt(x(i)^2+y(j)^2))^6-5/256*beta*(sqrt(x(i)^2+y(j)^2))^8));
         
         p=exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(j)^2));
         
         m=m+p;
    end
    test1=test1+m*dx*dy;
    
end
w=test1 % 计算一维p密度的主体积分的值；1/w ： 规范性常数

[X,Y]=meshgrid(x,y);
Pz=(1.0/w)*exp((-eta/sig2)*(alpha*X.^2+beta/2.0*X.^4+Y.^2));
%-----------------------------------------------------------------------
figure (1) 
mesh(X,Y,Pz);
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('联合概率密度f(x,y)')
view(58,60)

figure(2)
contour(x,y,Pz,70,'k')
xlabel('x')
ylabel('y')

%----------------------------------------------------------
for i=1:length(x)
    mid1=0;
    for j=1:length(y)
         p=(1.0/w)*exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(j)^2));
        mid1=mid1+p*dy;
    end
    px(i)=mid1;
end

testx=0;
for i=1:length(px)
    testx=testx+px(i)*dx;
end
testx
%----------------------------------------------------------
for j=1:length(y)
    mid1=0;
    for i=1:length(x)
         p=(1.0/w)*exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(j)^2));
        mid1=mid1+p*dx;
    end
    py(j)=mid1;
end

testy=0;
for j=1:length(py)
    testy=testy+py(j)*dy;
end
testy
%----------------------------------------------------------------
text_xy=0; k=100;   %计算给定y(k)处，x的条件概率密度函数
for i=1:length(x)   
    
        fk1(i)=(1.0/w)*exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(k)^2));
       fk1(i)=fk1(i)/py(k);
       text_xy=text_xy+fk1(i)*dx;
end
text_xy
y(k)
%----------------------------------------------------------------
text_xy=0; k=150;   %计算给定y(k)处，x的条件概率密度函数
for i=1:length(x)   
    
        fk2(i)=(1.0/w)*exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(k)^2));
       fk2(i)=fk2(i)/py(k);
       text_xy=text_xy+fk2(i)*dx;
end
text_xy
y(k)
%-------------------------------------------------
text_xy=0; k=350;   %计算给定y(k)处，x的条件概率密度函数
for i=1:length(x)       
       fk3(i)=(1.0/w)*exp((-eta/sig2)*(alpha*x(i)^2+beta/2.0*x(i)^4+y(k)^2));
       fk3(i)=fk3(i)/py(k);
       text_xy=text_xy+fk3(i)*dx;
end
text_xy
y(k)
%-------------------------------------------------


figure(3)
plot(x,px,'b-')
hold on 
% plot(x,fk1,'m--')
% hold on 
% plot(x,fk2,'k-.')
% hold on
% plot(x,fk3,'r-.')
% hold on
 xlabel('x')
 ylabel('f_X(x)')
% ylabel('f_X(x) 或 f_X_|_Y(x|y)')

