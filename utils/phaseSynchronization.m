function [h,c] = phaseSynchronization(ref, g)
% h = phaseSynchronization(ref, g)
% model: ref = c * g
% minimization of err = integral (f-c*g).^2 dx yields...

c = sum(sum( conj(g) .* ref )) / sum(sum( conj(g).*g ));
c = c / abs(c);

% h is phase synchronized
h = c * g;
end