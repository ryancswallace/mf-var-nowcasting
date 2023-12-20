function output = nancumsum(input)

ind    = ~isnan(input);
output = nan*input;
output(ind==1,:) = cumsum(input(ind==1,:));
