function WSample=Uiniform_Sampling(Poly,Num_sample)
% generate w_k samples
nx=size(Poly.A,2);
WSample=[];
Dmin=min(Poly.V,[],1);
Dmin=Dmin';
Dmax=max(Poly.V,[],1);
Dmax=Dmax';

for i=1:Num_sample
    flag=1;
    while flag
        wsample=(Dmax-Dmin).*rand(nx,1)+Dmin;
        if Poly.A*wsample<=Poly.b
            WSample=[WSample wsample];
            flag=0;
        end
    end
end

end