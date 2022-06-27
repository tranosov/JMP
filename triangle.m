function triangle=triangle(x1,x2,y1,y2,above,decr)
if (decr==1 && above==1) || (decr==0 && above==0) % incr and bellow also gives bigger weight to bigger values, is equivalent
    %lower=(eps1*(y2-y1)/(x1-x2) +y2 -x1*(y2-y1)/(x1-x2));
    if y2>y1 && x2>x1
    triangle=(2/((y2-y1)*(x2-x1)))*((y2-y1)/(x1-x2))*...
        (0.5*x1*(x2^2-x1^2) - ...
        (1/3)*(x2^3-x1^3) );
    else
        triangle=0;
    end
        
end
%todo: do the rest later
end