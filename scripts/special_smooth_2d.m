function [output] = special_smooth_2d(input, win, bins1, bins2, nGrid)

output = zeros(size(input));

nGrid1 = nGrid(1);
nGrid2 = nGrid(2);

x1 = (-nGrid1:nGrid1)/nGrid1;
x2 = (-nGrid2:nGrid2)/nGrid2;

Smoother1 = exp(-x1.^2/win(1)^2/2);
Smoother2 = exp(-x2.^2/win(2)^2/2);
Smoother1 = Smoother1./sum(Smoother1);
Smoother2 = Smoother2./sum(Smoother2);

delta1 = zeros(size(Smoother1)); delta1(nGrid1 + 1) = 1;
delta2 = zeros(size(Smoother2)); delta2(nGrid2 + 1) = 1;

if win(1)==1
    trash = nanmean(input,1);
    input1 = repmat(trash,[size(input,1) 1]);
elseif isnan(win(1))
    input1 = input;
else
    input1  = conv2(Smoother1, delta2,input, 'same');
end

if win(2)==1
    trash = nanmean(input1,2);
    output = repmat(trash,[1 size(input1,2)]);
elseif isnan(win(2))
    output = input1;
else
    output   = conv2(delta1,    Smoother2,  input1,'same');
end


end