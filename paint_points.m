A=importdata("2.62.txt");
x=A(:,1);y=A(:,2);z=A(:,3);
[X,Y]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y));
Z=griddata(x,y,z,X,Y,'v4');
surf(X,Y,Z)
shading interp
hold on;
plot(x,y,'.')
%plot3(x,y,z,'.');grid on;