function lp = logposterior(theta,X,Y)
                                                                            % theta(2) = alpha
                                                                            % theta(3) = m

% log_prior = log(normpdf(theta(2),12,4)) + log(normpdf(theta(3),0,10/365)); % log_prior = log(prod(pdf(par1)*pdf(par2)*...*pdf(parM)) = log(pdf(par1)) + log(pdf(par2)) + ...)
                                                                           % pdf(e0) = cost => log(pdf(e0)) neglected
                                                                           % pdf(sigma_LH) = cost => log(pdf(sigma_LH)) neglected
log_prior = log(normpdf(theta(2),12,4)) + log(normpdf(theta(3),0,10/365));


Yhat = mymodel(X, theta); % results of prediction model
log_lh = sum(log(normpdf(Y,Yhat,theta(4)))); %sum(log(pdf(Y, Yhat, sigma_LH));, Y = epsilon, Yhat = e0 + alpha*dT + m*d t
                                             % sum(log(pdf(z, 0, sigma_LH))
                                             % is the same. z = Y - Yhat
lp = log_prior + log_lh; % log_posterior = log_prior + log_likelihood  

end

