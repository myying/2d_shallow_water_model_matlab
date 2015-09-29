%solve Au=q for x with Thomas algorithm, given that A is tridiagonal.
function u=tridiagonal(A,q,n)
	for i=2:n
		r=A(i,i-1);
		A(i,:)=A(i,:)-r*A(i-1,:);
		q(i)=q(i)-r*q(i-1);
		r=A(i,i);
		A(i,:)=A(i,:)/r;
		q(i)=q(i)/r;
	end
	for i=n-1:-1:1
		r=A(i,i+1);
		A(i,:)=A(i,:)-r*A(i+1,:);
		q(i)=q(i)-r*q(i+1);
	end
	u=q;
end
