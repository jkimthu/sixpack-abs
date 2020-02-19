function s = sumRho(rho,t,f)

% given an array of T and f(T) values,
% 1. calculate the sum of all T and f(t) pairs, for each value of rho
% 2. output vector of sums

% initialize vector of output = sums as function of rho
sumsPerRho = zeros(length(rho),1);

for r = 1:length(rho)
    
    % initialize vector of solutions for each T and f(T) pair
    pairValues = zeros(length(t),1); 
    
    for i = 1:length(t)
        
        % calculate exponential for each T and f(T) pair
        pairValues(i) = exp(-rho(r)*t(i))*f(i);   
        
    end
    
    % calculate sum of solutions per rho
    sumsPerRho(r) = sum(pairValues);
    
end

% output all sums (row corresponds to that of input, rho)
s = sumsPerRho;

end

   
