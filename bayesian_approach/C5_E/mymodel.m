function eps = mymodel(X,THETA)

eps=ones(size(X,1),1)*THETA(1,1)+THETA(2,1)*X(:,1)+THETA(3,1)*X(:,2);

end

