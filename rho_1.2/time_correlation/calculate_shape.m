function [p,coord] = calculate_shape(rg_all)
%%
for i1 = 1:size(rg_all,3)
    [V,D] = eig(rg_all(:,:,i1));
    [lambda,ind] = sort(diag(D));
%     Vs(:,:,i1) = V(:,ind);
    b = lambda(3) - (lambda(1)+lambda(2))/2;
    c = (lambda(2)-lambda(1));
    asph(i1,:) = b;
    acyl(i1,:) = c;
    kappa(i1,:) = (b^2+c^2*3/4)/(lambda(1)+lambda(2)+lambda(3))^2;
    S6(i1,:) = (lambda(3)-lambda(1))/sqrt(2)/sum(lambda);
    S7(i1,:) = (2*lambda(2)-lambda(3)-lambda(1))/sqrt(6)/sum(lambda);
%     S6(i1,:) = (lambda(3)-lambda(1))/sqrt(2);
%     S7(i1,:) = (2*lambda(2)-lambda(3)-lambda(1))/sqrt(6);
end
p = [asph acyl kappa];
coord = [S6 S7];
end

