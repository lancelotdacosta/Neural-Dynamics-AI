function length = cat_inf_length(MDP)
length =0;
if(class(MDP) == 'struct')
    for agent = 1:size(MDP,1)
        for trial = 1:size(MDP,2)
            for state = 1:size(MDP(agent,trial).xn,2)
                length= length + ...
                    sum(sqrt(sum((sqrt(MDP(agent,trial).xn{1,state}(1:(end-1),:))...
                    -sqrt(MDP(agent,trial).xn{1,state}(2:end,:))).^2,2)));
            end
        end
    end
elseif(class(MDP) == 'double')
    length = sum(sqrt(sum((sqrt(MDP(1:(end-1),:))...
                    -sqrt(MDP(2:end,:))).^2,2)));
end
length= length*2;
    
        


function hello
x=[0.3,0.7;
    0.8,0.2];
cat_inf_leng(x)
fun = @(t) abs(1./sqrt(t.*(1-t))); 
integral(fun, x(1,2), x(2,2))
helloworld(x)
helloxxx(x)
helloy(x)

function hello = helloworld(x)
hello = sqrt(8*(1-sqrt(x(1,:))*sqrt(x(2,:))'));

function hello = helloxxx(x)
hello = 2*acos(sqrt(x(1,:))*sqrt(x(2,:))');

function hello = helloy(x)
hello = 2*norm(sqrt(x(1,:))-sqrt(x(2,:)));



function length = cat_inf_leng(y) 
%loop over size of x to get list of updates

%makes sure that two consecutive updates which are the same are not
%computed
% ConsecDuplic  = [1 2 2 3 4 5 5 5 6 7 8 9 9]
% LV = [false ConsecDuplic(2:end) == ConsecDuplic(1:end-1)]
% NoDuplic = ConsecDuplic(~LV)

%only if remaining array has at least two rows
x=y;
%from now on suppose x is my list of updates
a = x(1:end-1,:)-2*sqrt(x(1:end-1,:).*x(2:end,:))+x(2:end,:);
length = 2*sum(sqrt(sum(a,2)));



% function length = cat_inf_lengt(y) %got a problem as fun returns a tensor
% %loop over size of x to get list of updates
% 
% %makes sure that two consecutive updates which are the same are not
% %computed
% % ConsecDuplic  = [1 2 2 3 4 5 5 5 6 7 8 9 9]
% % LV = [false ConsecDuplic(2:end) == ConsecDuplic(1:end-1)]
% % NoDuplic = ConsecDuplic(~LV)
% 
% %only if remaining array has at least two rows
% x=y;
% %from now on suppose x is my list of updates
% a = x(1:end-1,:)-2*sqrt(x(1:end-1,:).*x(2:end,:))+x(2:end,:);
% b = 2*(-x(1:(end-1),:)+sqrt(x(1:end-1,:).*x(2:end,:)));
% gamma0 = @(t) spm_cross(a,t.^2)+spm_cross(b,t)+spm_cross(x(1:end-1,:),ones(1,numel(t))); %geodesic path % t are line vectors, a is column vector
% gamma1 = @(t) spm_cross(2*a,t)+spm_cross(b,ones(1,numel(t))); %derivative of geodesic
% fun = @(t) sum(reshape(sqrt(sum(gamma1(t).^2./gamma0(t),2)), [size(a,1) 1 numel(t)]),1); 
% t=(1:9)/10;
% length = sqrt(integral(fun,0,1));  %information length

function length = categorical_inf_length(d) 
x=d(1,:)';
y=d(2,:)';
a = x-2*sqrt(x.*y)+y;
b = 2*(-x+sqrt(x.*y));
gamma0 = @(t) a*t.^2+b*t+x*ones(1,numel(t)); %geodesic path % t are line vectors, a is column vector
gamma1 = @(t) 2*a*t+b*ones(1,numel(t)); %derivative of geodesic
fun = @(t) sum(gamma1(t).^2./gamma0(t)); 
%fun((1:9)/10)
length = sqrt(integral(fun,0,1));