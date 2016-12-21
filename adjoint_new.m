function [y, dy] = adjoint_new(k, x_0, x)
C0 = (exp(1/k)-exp(x_0/k))/(1-exp(1/k));
C1 = -C0;
C2 = (exp(1/k)-exp((1+x_0)/k))/(1-exp(1/k));
C3 = -((1-exp(-x_0/k))/(1-exp(1/k)));

for i = 1:numel(x)
	if (x(i) <= x_0)
		K1 = C0;
		K2 = C1;
	else
		K1 = C2;
		K2 = C3;
	end

	y(i) = K1 * exp(-x(i)/k) + K2;
	dy(i) = - K1 * exp(-x(i)/k) / k;
end
end
