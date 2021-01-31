%b(1)=b^0_k,b(2)=b^1_k
%Need b(2)>=b(1)
%c(1)=c_1+,c(2)=c_1-, c(3)=c_2
%y1 is the equation of motion using c_1+, y2 uses c_2-

%Set up parameters
b(1)=0.1;
b(2)=0.6;
c(1)=-2*b(1)*(1+sqrt(1+(b(2)-b(1))/b(1)));
c(2)=-2*b(1)*(1-sqrt(1+(b(2)-b(1))/b(1)));
c(3)=b(1);

%Equation of motion

y1 = @(t) (c(1)*t).^2/(4*c(3))+c(1)*t+c(3);
y2 = @(t) (c(2)*t).^2/(4*c(3))+c(2)*t+c(3);

n =10;
t= (0:n)/n;

%plot(t,y1(t),'Color',[0 1 0]), hold on
%plot(t,y2(t),'Color',[1 0 0]), hold off %works
a = y2(t);

%Set up parameters
b(1)=1-b(2);
b(2)=1-b(1);
c(1)=-2*b(1)*(1+sqrt(1+(b(2)-b(1))/b(1)));
c(2)=-2*b(1)*(1-sqrt(1+(b(2)-b(1))/b(1)));
c(3)=b(1);

y2 = @(t) (c(2)*t).^2/(4*c(3))+c(2)*t+c(3);

b = y2(t);

a
b
