

val = 0;
xgrad = grad;
for j = 1:ndat
    if x(j) <= glaveps
        if xgrad(j) < 0 
            val = max(val, -xgrad(j)); 
        end
%     elseif x(j) >= C(j)%-glaveps
%         if xgrad(j) > 0 val = max(val, xgrad(j)); end
    else
        val = max(val,abs(xgrad(j)));
    end
end
%val = norm(grad,inf);