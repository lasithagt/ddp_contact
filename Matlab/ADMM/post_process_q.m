q_ = q;
for i = 1:size(q,1)
    mask = q(i,:)>pi;
    q_(i,mask) =  q(i,mask) - 2*pi;
end