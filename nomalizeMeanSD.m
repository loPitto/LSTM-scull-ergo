function [OUT, mu, sigma] = nomalizeMeanSD(IN)
mu    = mean(cat(2,IN{:}),2);
sigma = std(cat(2,IN{:}),0,2);
for n = 1:numel(IN)
    OUT{n} = (IN{n} - mu) ./ sigma;
end