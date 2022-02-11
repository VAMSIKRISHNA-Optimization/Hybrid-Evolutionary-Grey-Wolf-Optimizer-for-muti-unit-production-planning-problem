

%%Hybrid Evolutionary Grey Wolf Optimizer for muti-unit production planning problem
% By Aala Kalananda Vamsi Krishna Reddy 

function [bestfitness,bestsol,TIM] =HEGWO(Np,T,lb,ub,dim,fobj,func,Fact)
tic;
%% Starting of HE-GWO
%Fact is the fitness value of the global best            
f = NaN(Np,1);             % Vector to store the fitness function value of the population members
fu = NaN(Np,1);          % Vector to store the fitness function value of the New population members
D = dim;                        % Determining the number of decision variables in the problem
U = NaN(Np,D);          % Matrix to store the trial solutions




P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);   % Generation of the initial population
for p = 1:Np
%     f(p) = fobj(P(p,:));           % Evaluating the fitness function of the initial population
    f(p)=abs(Fact-feval(fobj,P(p,:)',func)); %cec2020
end

%% Iteration loop
for t = 1: T

     a=2-t*((2)/T)+0.05; % a decreases linearly fron 2 to 0 %%dependance factor

    for i = 1:Np
    
    
  

    Npcc=Gen_cr(i,Np);  %%Compliance rates


        %% Hunting
        Candidates = [1:i-1 i+1:Np];            % Ensuring that the current member is not the partner
        rdx = Candidates(randperm(Np-1,6));     % Selection of three random parters
        
        [~,idx]=sort(f);

        X1 = P(idx(1),:);                %%Alpha Wolf        
        X2 = P(idx(2),:);                       %%Beta Wolf
        X3 = P(idx(3),:);                               %%Delta Wolf
        

        
        
        R1 = P(rdx(1),:);                       % Assigning randomly selected solution 1
        R2 = P(rdx(2),:);                       % Assigning randomly selected solution 2
        R3 = P(rdx(3),:);                       % Assigning randomly selected solution 3
          R4 = P(rdx(4),:);                      % Assigning randomly selected solution 4
      
        



        
        V1 =  X1 + a*(R1 - R2);
            V2 =  X2 + a*(R3 - R4); 
                    V3 =  X3;

        V=(V1+V2+V3)/3;
       
        
        %% Crossover




        del = randi(D,1);                        % Generating the random varibale delta
        for j = 1:D

            if (rand <= Npcc) || del == j          % Check for donor vector or target vector
                U(i,j) = V(j);                   % Accept variable from donor vector
            else
                U(i,j) = P(i,j);                 % Accept variable from target vector
            end
            
        end
        
        


%     %% Bounding and Greedy Selection

        
        U(i,:) = min(ub,U(i,:));                % Bounding the violating variables to their upper bound
        U(i,:) = max(lb,U(i,:));                % Bounding the violating variables to their lower bound
        

        fu(i)=abs(Fact-feval(fobj,U(i,:)',func)); %cec2020

        if fu(i) < f(i)                         % Greedy selection
            P(i,:) = U(i,:);                    % Include the new solution in population
            f(i) = fu(i);                       % Include the fitness function value of the new solution in population;
            

         
        end
    end
    


%%SHUFFLING THE HUNTERS

     sdx = Candidates(randperm(Np-1,Np-4));

for i=1:size(Np,1)
        P(i,:)=P(sdx(i),:);
end

    
    
    [bestfitness,ind] = min(f);
    TIM=toc;

end
[bestfitness,ind] = min(f);
bestsol = P(ind,:);


end

%%Function to generate compliance rates
function [Npcc] = Gen_cr(i,Np)

if i<round(Np*0.25)
    Npcc=(0.2-0.1).*rand(1,1)+0.1;
elseif i>=round(Np*0.25) &&i<=round(Np*0.5)
    Npcc=(0.4-0.1).*rand(1,1)+0.1;
    elseif i>=round(Np*0.5) &&i<=round(Np*0.75)
    Npcc=(0.6-0.1).*rand(1,1)+0.1;
else
    Npcc=(0.8-0.1).*rand(1,1)+0.1;
end

