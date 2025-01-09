function answer = trace_sum_nabla(u, L)
answer = zeros(size(L,2),size(L,3));
for i = 1:size(L,3)
    answer(:,i) = L(:,:,i)*u(:,i);
end
answer = answer*2;