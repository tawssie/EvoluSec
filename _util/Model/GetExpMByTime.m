function M=GetExpMByTime(V,D,Time)
% ExpM means exponential matrix
% This function generates exponential matrices 
% from eigenvectors and eigenvalues of the evolutionary model.
%
% The exponential matrix therefore represents transistion probability at
% "Time" according to the evolutionary model.
%
% V = matix of eigenvectors
% D = 1. square matrix with eigenvalue on diagonal
%     2. can be a vector of eigenvalues
%
% if Time is a single value, M is the square probability matrix.
% if Time is an array, then M is a cell list with all corresponding probability matrix

Dim=size(D);
if Dim(1)==Dim(2)&&Dim(1)~=0
    EigenValue=diag(D);
else
    EigenValue=D;
end


if numel(Time)>1
    
    MList=cell(numel(Time),1);
    
    for i=1:numel(Time)
        
        if Time(i)==0 %%% to prevent precision error by assigning identity matrix directly
            M=eye(size(V));
        else
        
            M=V*diag(exp(EigenValue*Time(i)))/V;

            % Warning any negative probability because of precesion
            if any(M(:)<0)
                warning('prob < 0 !? at the time %.10f',Time(i));
            end
        end

        MList{i}=M;
    end
    M=MList;
else

    if Time==0 % Prevent precision issues by assigning identity matrix directly at zero Time.
        M=eye(size(V));
    else
    
        M=V*diag(exp(EigenValue*Time))/V;
        % Warning any negative probability because of precesion
        if any(M(:)<0)
            warning('prob < 0 at the time %.10f',Time);
        end
    end

end
