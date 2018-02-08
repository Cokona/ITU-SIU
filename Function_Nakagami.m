% Step size sabit (epsilon)
% Rayleigh

function [Probability_Container,Average_It_Container] = Function_Nakagami(SNR)

%% Defining the System
MaxIt = 2000;
epsilon = 1/100;   % step size determined random
min_node = 4;
max_node = 7;
Allowed_Error = 1/1000; 
nodes = min_node:1:max_node;

for n = 1:length(nodes)
    
    [ii,jj] = ndgrid(1:nodes(n));              % Used to choose the upper triangle of The Matrix
    Edge_num = nchoosek(nodes(n),2);           % for us to be able to define any number of edges
    monte_counter = 0;
    counter = 0;
    switch nodes(n)
%       case 3
%           Measured = [-6 ; -3 ; 7];   
      case 4
          Measured = [-6 ; -3 ; 7 ; 14 ];
      case 5
          Measured = [-6 ; -3 ; 7 ; 14; 21];
      case 6
          Measured = [-6; -3; 7; 14; 21; -12.5];
      case 7
          Measured = [-6; -3; 7; 14; 17; -12.5; 8.5];
%       case 8
%           Measured = [-6; -3; 7; 14; 17; -12.5; 8.5; -9.5];
    end     
IterMeasured = Measured;

% montemax = 1000 ;
    if SNR < 3
    montemax = 5e5; 
    elseif SNR < 6
    montemax = 1e6;
    elseif SNR < 10
    montemax = 5e6;
    elseif SNR < 15
    montemax = 1e7;
    else
    montemax = 1e7;
    end
  
    Ranks = zeros(montemax,3);                %this will keep everything about matrices
    
    %% Starting to Monte Carlo
    
    for monte = 1:montemax
     %Nakagami baþlangýcý
        R = 1 ; 
        Naka_m=2;
        M=100000;
        h=[];
    
        for d=1: 2*Naka_m
            h(d,:)= (((randn(1,Edge_num))/sqrt(2*Naka_m)));

        end
        h=h.^2;
      
        if Naka_m==0.5    %%Gauss
            h_Naka=sqrt(h);
        else
            h_Naka=sqrt( sum(h));
        end

        gama = 10^(SNR/10); 
        threshold = (2^R - 1) / gama ;
        Edge_con = abs(h_Naka).^2 > threshold ;

      A = zeros(nodes(n));                       % Adjacency Matrix's for all
      A(jj>ii) =  Edge_con;                   % Filling the upper Triangle
      %WOW! A good way of filling a upper triangular according to threshold
      
      A = A + A';                             % Adjacency Matrix created
      D = diag(sum(A));                       % In-Degree Matrix created  
      L = D - A ;                             % Laplacian Matrix
     
      Ranks(monte,1) = rank(L) ;  %her bir monte de?i?keni için matrix rank toplamaca
      
      %% Now let's start iterating for every matrix in monte carlo
      
      if  Ranks(monte,1) == nodes(n) - 1 %bu sa?lanm?yorsa alphaya ula??lamad? ranks5 = 0
     
      Ranks(monte,2) = 1; %we can reach alpha,consensus
      monte_counter = monte_counter + 1;
      
      for k = 1:MaxIt

      P_epsilon = eye(nodes(n)) - epsilon * L;
      IterMeasured = P_epsilon * IterMeasured;  
      
      %% Defining an Error Rule
      if abs(max(IterMeasured) - min(IterMeasured)) > Allowed_Error 
      counter = counter + 1;
      else
          Ranks(monte,3) = counter;
          IterMeasured = Measured;
          counter = 0;
          break 
      end
      
      end
      
      else
      Ranks(monte,2) = 0;
      Ranks(monte,3) = 0;
      end
      
end
%% Calculating the Probabilities of Ranks

   Average_It_per_SNR = sum(Ranks(:,3))/sum(Ranks(:,2));
   Average_It_Container(1,n) = Average_It_per_SNR;
   Succ_Probability = sum(Ranks(:,2)) / montemax;
   Probability_Container(1,n) =  Succ_Probability;

end