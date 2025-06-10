function SNodeLeafA=SetTransitionProb(isPrior,SNodeLeafA,Domain,MA0)
% isPrior == true  => assign prob as prior (stationary prob/normalised frequencies)
% isPrior == false => assign prob with transition matrix (expm)
import bn.prob.EnumDistrib;

if isPrior
    SNodeLeafA.put(EnumDistrib(Domain, MA0));
else
    AAOrder=Domain.getValues; % A,C,D...Y

    for i=1:numel(AAOrder)
        SNodeLeafA.put(EnumDistrib(Domain, MA0(i,:)),AAOrder(i));
    end
end
