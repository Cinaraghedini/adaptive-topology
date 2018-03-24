
close all

figure(18)

results=zeros(20,5);

for i=1:3 
    pd_gamma = makedist('Gamma','a',4,'b',4)
    pd_gamma=truncate(pd_gamma,1,80);
    gamma_rand_numbers = random(pd_gamma,20,1);
    results(:,1)=results(:,1)+gamma_rand_numbers;
    pd_gamma = makedist('Gamma','a',8,'b',4);
    pd_gamma=truncate(pd_gamma,1,80);
    gamma_rand_numbers = random(pd_gamma,20,1);
    results(:,2)=results(:,2)+gamma_rand_numbers;
    pd_gamma = makedist('Gamma','a',12,'b',4);
    pd_gamma=truncate(pd_gamma,1,80);
    gamma_rand_numbers = random(pd_gamma,20,1);
    results(:,3)=results(:,3)+gamma_rand_numbers;
    gamma_rand = sort(round([(1 - 80).*rand(20,1) + 80]));
    results(:,4)=results(:,4)+gamma_rand;
end
results=results/i;
figure(1)
[f,xi] = ksdensity(results(:,1));
plot(xi,f,'r');
hold all;
[f,xi] = ksdensity(results(:,2));
plot(xi,f,'b');
hold all;
[f,xi] = ksdensity(results(:,3));
plot(xi,f,'m');
hold all;
[f,xi] = ksdensity(results(:,4));
plot(xi,f,'g');
hold all;













pd_gamma = makedist('Gamma','a',4,'b',4);
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
results(:,1)=results(:,1)+gamma_rand_numbers;
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'r');
hold all;
pd_gamma = makedist('Gamma','a',8,'b',4)
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'b');
hold all;
pd_gamma = makedist('Gamma','a',12,'b',4)
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'m');
hold all;
gamma_rand = sort(round([(1 - 80).*rand(20,1) + 80]));
mean(gamma_rand)
std(gamma_rand)
[f,xi] = ksdensity(gamma_rand);
plot(xi,f,'g');


%%%

figure(3)

pd_gamma = makedist('Normal','mu',10,'sigma',10) 
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'r');
hold all;
pd_gamma = makedist('Normal','mu',30,'sigma',10) 
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'b');
hold all;
pd_gamma = makedist('Normal','mu',60,'sigma',10) 
mean(pd_gamma)
std(pd_gamma)
pd_gamma=truncate(pd_gamma,1,80);
gamma_rand_numbers = random(pd_gamma,20,1);
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'m');
hold all;
gamma_rand_numbers = sort(round([(1 - 80).*rand(20,1) + 80]));
mean(gamma_rand_numbers)
std(gamma_rand_numbers)
[f,xi] = ksdensity(gamma_rand_numbers);
plot(xi,f,'g');



