% Step size sabit (epsilon)
% Rayleigh

clear; clc; close all; %dbstop if error;
SNR = 0:20;

cores = 6;
c = parcluster('local');  % to start parallel pools
c.NumWorkers=cores ;
parpool(c, c.NumWorkers);


parfor s = 1:length(SNR)
    tic
    [Probability,Average_It] = Function_Rayleigh(SNR(s));
    
    Probability_Container(:,s) = Probability.';
    Average_It_Container(:,s) = Average_It.';
    toc
end

% Probability_Container = squeeze(Probability_Container);
% Average_It_Container = squeeze(Average_It_Container);

delete(gcp('nocreate'));   % to close parallel pools


number = [4:7]';
names = int2str(number);
                         
figure;
semilogy(SNR,1-Probability_Container(1,:),'x -');
hold on;
for k = 2:size(Probability_Container,1)
    semilogy(SNR,1-squeeze(Probability_Container(k,:)),'x -');
end
legend(names);
title('Propability of NOT Reaching Consensus');
savefig('Rayleigh_fig1.fig')
hold off

figure;
plot(SNR,Average_It_Container(1,:),'x -');
hold on;
for k = 2:size(Probability_Container,1)
    plot(SNR,Average_It_Container(k,:),'x -');
end
legend(names);
title('Average Number of Iterations Until Consensus Is Reached');
savefig('Rayleigh_fig2.fig')


clear c ans
save('Rayleigh_data_with0.mat')

