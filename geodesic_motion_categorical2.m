%Set up parameters
x(1)=0.001;
x(2)=0.999;


%Equation of motion
e = 1;
y1 = @(t) (x(1)-2*e*sqrt(x(1)*x(2))+x(2))*t.^2+2*(-x(1)+e*sqrt(x(1)*x(2)))*t+x(1);

e = -1;
y2 = @(t) (x(1)-2*e*sqrt(x(1)*x(2))+x(2))*t.^2+2*(-x(1)+e*sqrt(x(1)*x(2)))*t+x(1);

n =1000;
t= (0:n)/n;

plot(t,y1(t),'Color',[0 1 0]), hold on %works
plot(t,y2(t),'Color',[1 0 0]), hold off
%implies epsilon =1

