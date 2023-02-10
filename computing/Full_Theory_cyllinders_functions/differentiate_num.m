function [df] = differentiate_num(f,x)
% differentiate_num numerically evaluates the derivative of a univariate 
% function.

N = size(x,2)-1;

% preallocates derivative results
df = zeros(size(f), 'gpuArray');

% derivative at lower bound using forward difference
df(:,1,:) = (f(:,2,:)-f(:,1,:))./(x(:,2,:)-x(:,1,:));

% derivative at upper bound using backward difference
df(:,N+1,:) = (f(:,N+1,:)-f(:,N,:))./(x(:,N+1,:)-x(:,N,:));

% derivatives at all other nodes using central differences
df(:,2:N,:) = (f(:,3:end,:)-f(:,1:end-2,:))./(x(:,3:end,:)-x(:,1:end-2,:));

end