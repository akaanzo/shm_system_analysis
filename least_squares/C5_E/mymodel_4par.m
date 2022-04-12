function eps=mymodel_4par(X,theta)

%eps=ones(size(t))*theta(1)+theta(2)*t+theta(3)*T;
%X=[t T];
%theta=[epsilon0; m; alpha]; 
eps=ones(size(X,1),1)*theta(1)+theta(2)*X(:,1)+theta(3)*X(:,2)+theta(4)*X(:,2);

end

