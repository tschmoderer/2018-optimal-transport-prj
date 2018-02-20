function dctf = mydct2(f,N,Q)
dctf = zeros(Q+1,N+1);
t    = zeros(Q+1,N+1);
for i = 1:Q+1
   t(i,:) = mydct(f(i,:),N);
end
for i = 1:N+1
  dctf(:,i) = mydct(t(:,i)',Q)';  
end

%dctf = mydct(mydct(f,


%b = dct(dct(arg1).').';
%     for i = 1:Q+1
%         for j = 1:N+1
%             s = 0;
%             for k = 1:Q+1
%                 for l = 1:N+1
%                     s = s + f(k,j)*cos(pi*(2*k-1)*(i-1)/(2.*(Q+1)))*cos(pi*(2*l-1)*(j-1)/(2.*(N+1)));
%                 end 
%             end
%             if (i == 1) 
%                 ai = 1/sqrt(Q+1);
%             else 
%                 ai = sqrt(2/(Q+1));
%             end
%              if (j == 1) 
%                 aj = 1/sqrt(N+1);
%             else 
%                 aj = sqrt(2/(N+1));
%              end
%             dctf(i,j) = ai*aj*s;
%         end
%     end
end

