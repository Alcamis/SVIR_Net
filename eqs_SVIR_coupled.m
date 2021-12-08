function [Y]=eqs_SVIR_coupled(a1,a2,N,t,x,C,L)

	mu=C(:,1);
    beta=C(:,2);
    phi=C(:,3);
    rho=C(:,4);
    lambda=C(:,5);
    delta=C(:,6);
    theta=C(:,7);

    % S: x(1+4*(j-1))
    % V: x(2+4*(j-1))
    % I: x(3+4*(j-1))
    % R: x(4+4*(j-1))
    
    for j=1:N % loop over all N communities
        ss=0;
        si=0;
        for m=1:N % loop over all N communities
            ss=ss+L(j,m)*x(1+4*(m-1));
            si=si+L(j,m)*x(3+4*(m-1));
        end
        Y(1+4*(j-1))=mu(j)-beta(j).*x(1+4*(j-1)).*x(3+4*(j-1))-phi(j).*x(1+4*(j-1))-mu(j).*x(1+4*(j-1))+theta(j).*x(2+4*(j-1))+delta(j).*x(4+4*(j-1))-a1.*ss;
        Y(2+4*(j-1))=phi(j).*x(1+4*(j-1))-rho(j).*beta(j).*x(2+4*(j-1)).*x(3+4*(j-1))-mu(j).*x(2+4*(j-1))-theta(j).*x(2+4*(j-1));
        Y(3+4*(j-1))=beta(j).*x(1+4*(j-1)).*x(3+4*(j-1))+rho(j).*beta(j)*x(2+4*(j-1)).*x(3+4*(j-1))-lambda(j).*x(3+4*(j-1))-mu(j).*x(3+4*(j-1))-a2.*si;
        Y(4+4*(j-1))=lambda(j).*x(3+4*(j-1))-mu(j).*x(4+4*(j-1))-delta(j).*x(4+4*(j-1));
    end
    Y=Y';
end