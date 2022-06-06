function concOut = diffTimestepCN(concIn,dt,dx,D,periodic)
%DIFFTIMESTEPCN applies the Crank-Nicholson method to simulate 2D diffusion,
%using the fractional steps implementation. See Cen et al. (2016)
%
%   INPUTS:
%       -concIn: NxM matrix containing the current values of the
%       concentration field
%       -dt: The timestep size
%       -dx: The pixel separation
%       -D: The diffusion constant
%       -periodic: Whether periodic boundary conditions should be applied
%
%   OUTPUTS:
%       -concOut: The updated concentration field
%
%   Author: Oliver J. Meacock (c) 2021

A = (dt*D)/(2*(dx^2));

%% The x-sweep
halfSol = zeros(size(concIn));

%Solve RHS
xRHS = (1-2*A)*concIn + A*circshift(concIn,1,2) + A*circshift(concIn,-1,2);

if ~periodic
    xRHS(:,1) = xRHS(:,1) - A*concIn(:,end);
    xRHS(:,end) = xRHS(:,end) - A*concIn(:,1);
end

%Solve tridiagonal matrix constructed from each row of the input
%concentrations and the RHS
if ~periodic
    a = ones(size(concIn,2),1) * -A;
    b = ones(size(concIn,2),1) * (1 + 2*A);
    c = a;
    a(1) = 0;
    c(end) = 0;
else
    solMat = eye(size(concIn,2));
    solMat = solMat * (1 + 2*A) - A*(circshift(solMat,1,1) + circshift(solMat,-1,1));
    solMatInv = inv(solMat);
end

for i = 1:size(concIn,1)
    %If non-periodic BCs, can use a much more computationally cheap
    %function to calculate necessary matrix
    if ~periodic
        halfSol(i,:) = tridiagSparseSolve(a,b,c,xRHS(i,:)');
    else
        halfSol(i,:) = solMatInv*xRHS(i,:)';
    end
end

%% The y-sweep
concOut = zeros(size(concIn));

%Solve RHS
yRHS = (1-2*A)*halfSol + A*circshift(halfSol,1,1) + A*circshift(halfSol,-1,1);

if ~periodic
    yRHS(1,:) = yRHS(1,:) - A*halfSol(end,:);
    yRHS(end,:) = yRHS(end,:) - A*halfSol(1,:);
end

%Solve tridiagonal matrix constructed from each row of the input
%concentrations and the RHS
if ~periodic
    a = ones(size(concIn,1),1) * -A;
    b = ones(size(concIn,1),1) * (1 + 2*A);
    c = a;
    a(1) = 0;
    c(end) = 0;
else
    solMat = eye(size(concIn,1));
    solMat = solMat * (1 + 2*A) - A*(circshift(solMat,1,1) + circshift(solMat,-1,1));
    solMatInv = inv(solMat);
end

for i = 1:size(concIn,2)
    if ~periodic
        concOut(:,i) = tridiagSparseSolve(a,b,c,yRHS(:,i));
    else
        concOut(:,i) = solMatInv*yRHS(:,i);
    end
end