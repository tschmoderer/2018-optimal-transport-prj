function f = myidct2(df,N,Q)
    f = zeros(Q+1,N+1);
    t = zeros(Q+1,N+1);
    for i = 1:Q+1
       t(i,:) = myidct(df(i,:),N);
    end
    for i = 1:N+1
      f(:,i) = myidct(t(:,i)',Q)';  
    end
end

