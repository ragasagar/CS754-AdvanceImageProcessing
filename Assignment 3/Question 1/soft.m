function theta = soft(y, lambda)
    theta = zeros(size(y));
    for i=1:length(y)
        if y(i) >= lambda
            theta(i) = y(i) - lambda;
        elseif y(i)<= -lambda
            theta(i) = y(i)+lambda;
        else
            theta(i) = 0;
        end
    end