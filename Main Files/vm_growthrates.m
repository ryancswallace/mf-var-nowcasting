function output = vm_growthrates(input,select)

YY = exp(input);
YY(:,select==1) = input(:,select==1);

output = YY(2:end,:)./YY(1:end-1,:) - 1;

