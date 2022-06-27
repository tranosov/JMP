function Jc=logitmatching(Jw,Jm,Z)
I=size(Z);
I=I(1);
Jc=zeros(I,I);
if Jw(1)~=1 && Jw(2)~=1 && Jm(1)~=1 && Jm(2)~=1
    for i=1:I %men are is, women are js
        for j=1:I
            Jc(i,j)=exp(Z(i,j)/2-Z(1,1)/2 +log(Jw(j))+log(Jm(i))-log(Jw(1))-log(Jm(1)));
        end
    end
else
    Jc=Jm'*Jw;
end
Jc=Jc/sum(sum(Jc));



end