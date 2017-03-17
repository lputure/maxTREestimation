function UB = myBoundLBWMEC(B,P,Ns)
% sample Ns points (x,y,z) from half unit sphere (z>0) and choose the largest value

% input
% B: 3*Nf matrix, Nf fiducial points in 3D space
% P: 3*Np matrix, Np target locations in 3D space
% Ns: sample number of the rotation axis

% output
% UB: 1*Np array, the estimated maximum TRE magnitude of the Np target locations

[~,Np] = size(P);
SampleAxis = randn(3,Ns);
SampleAxis = SampleAxis./repmat(sqrt(sum(SampleAxis.*SampleAxis)),3,1); % normalized to unit length

%% [Y,V,F] = myLBWMEC(B,P)

Record = zeros(Ns,Np);
for ii = 1:Ns
    n = SampleAxis(:,ii); % take our an axis
    [U,~,~] = svd(n);
    Bx = U(:,2)'*B;
    By = U(:,3)'*B; 
    Px = U(:,2)'*P;
    Py = U(:,3)'*P; % Project to the orthogonal space
%     % log file
%     fid = fopen('log.txt','a');
%     fprintf(fid,'%d\n',ii);
%     fprintf(fid,'%f\n',[Bx;By]);
%     fprintf(fid,'%f\n',[Px;Py]);
%     fclose(fid);
    [~,V,~] = myLBWMEC([Bx;By],[Px;Py]); % solve the Localization Based Weighted Minimum Enclosing Circle Problem
    Record(ii,:) = V;
end
UB = max(Record);
