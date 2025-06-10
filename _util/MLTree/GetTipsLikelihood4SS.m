function [TipLikelihood,MSALabel]=GetTipsLikelihood4SS(MSSA,Model)

AcceptableChar=[Model.AAOrder(:)' '_-']; % BEGHIST_-

if any(~ismember(MSSA(:),AcceptableChar))
    error('unexpected character occurs, only accepts BEGHIST-_, terminated');
end


[~,MSALabel]=ismember(MSSA,AcceptableChar); % transform into index according to AA order

% MSALabel stands for different character following Model.AAOrder
% MSALabel == 1 <=> B
% MSALabel == 2 <=> E
% MSALabel == 3 <=> G
% MSALabel == 4 <=> H
% MSALabel == 5 <=> I
% MSALabel == 6 <=> S
% MSALabel == 7 <=> T
% MSALabel == 8 <=> _
% MSALabel == 9 <=> -

% TipLikelihood is 7*9 matrix, pre calculated to speed up
% the last column vector (MSALabel==9) stands for the gap.

TipLikelihood=[eye(7) Model.F(:) ones(7,1)]; 