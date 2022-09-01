N = 6;                      % Population size 
n = 2;                      % Number of variables
nb = 5;                     % Number of bits to represent each variable
Pc = 1;                     % Cross_over probability
Xmax = 0.5;                 % Geometric constraints
Xmin = 0;
Gmax = 10;                  % Maximum number of generations

t = zeros(N,n*nb);           % Matrix containing the solutions
for i=1:N                    
    for j=1:(n*nb)
        if( rand()<0.5 )
           t(i,j) = 0;
        else
           t(i,j) = 1;
        end
    end
end

Gn = 0;
x = zeros(n,1);
x_t = zeros(1,n);
avg_fit = zeros(1,Gmax);
max_fit = zeros(1,Gmax);
min_fit = zeros(1,Gmax);
min_x = zeros(Gmax,n);
max_x = zeros(Gmax,n);

while Gn < Gmax            
fit = zeros(N,1);           
for i=1:N                    
    m = zeros(n,1);
    for k=1:n
        for j=1:nb
            m(k,1) = m(k,1) + power(2,nb-j)*t(i,j+((k-1)*nb) );
        end
        x_t(1,k) = Xmin + ( (Xmax - Xmin)*m(k,1) )/( power(2,nb)-1 );
    end

    fit(i) = objfunc(x_t);
end

sum=0;          %Reproduction
for i=1:N
    sum = sum + fit(i);
end
K = 1/sum ;

P = zeros(N,1);         
for i=1:N
    P(i) = K*fit(i);
end

RW = zeros(N,1);            %Roulette wheel with angle of each element 
for i=1:N
   if(i==1)
       RW(i) = P(i)*360;
   else
       RW(i) = RW(i-1) + P(i)*360;
   end
end

pool = zeros(N,n*nb);      %Mating pool solutions

for i=1:N                  %Mating pool selection
    r = randi(360);
    for j=1:N
        if( r < RW(j) )
            index = j;
            break;
        end
    end
    pool(i,:) = s(j,:); 
end

par = zeros((N/2),3);    %Crossover   
temp = randperm(N,N);    
for i=1:N/2              
   par(i,1) = temp(i);
   par(i,2) = temp(i + (N/2));
   par(i,3) = rand();
end

rr = zeros(1,3);
rr(1,:) = randperm((n*nb-1),N/2);
for i=1:N/2                %Checking crossover probability
    if( par(i,3) <= Pc )
        r1 = rr(1,i);
        temp = zeros(1,(n*nb-r1));
        temp(1,:) = pool( par(i,1) , r1+1 : n*nb );          
        pool( par(i,1) , r1+1 : n*nb ) = pool( par(i,2) , r1+1 : n*nb );    %Cross over
        pool( par(i,2) , r1+1 : n*nb ) = temp(1,:);
    end
end

Gn = Gn +1;

avg = 0;                %Contains the average fitness avlue
max = -1000;
min = 1000;

for i=1:N                    
    m = zeros(n,1);
    for k=1:n
        for j=1:nb
            m(k) = m(k) + power(2,nb-j)*t(i,j+((k-1)*nb) );
        end
        x_t(1,k) = Xmin + ((Xmax - Xmin)*m(k))/( power(2,nb)-1 );
    end
    
    fitness = objfunc(x_t);
    
    if( fitness < min )
        min = fitness;
        min_x(Gn,:) = x_t;
    end
    
    if( fitness > max )
        max = fitness;
        max_x(Gn,:) = x_t;
    end
    
    avg = avg + fitness;
end

avg = avg /N ;

avg_fit(Gn) = avg;
max_fit(Gn) = max;
min_fit(Gn) = min;

end                         %% End of while loop

sol = zeros(N,n);           %% matrix haing the final solutions

for i=1:N                   %% loop to calculate final solution    
    m = zeros(n,1);
    for k=1:n
        for j=1:nb
            m(k) = m(k) + power(2,nb-j)*t(i,j+((k-1)*nb) );
        end
        x_t(1,k) = Xmin + ((Xmax - Xmin)*m(k))/( power(2,nb)-1 );
        sol(i,k) = x_t(1,k);
    end
end

fi = fopen('OUTPUT','w'); 
fprintf(fi,'\nThe value of ( x1, x2 ) having maximum fitness after each generation  %d\n');
fprintf(fi,'\n');
for i=1:Gmax
    for j=1:n
        fprintf(fi,'%d\t',max_x(i,j));
    end
    fprintf(fi,'\n');
end

    
    