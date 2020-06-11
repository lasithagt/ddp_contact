function p = convert_q2x(x)

    q = x(1:7,:);
    [Slist, M] = manipulator_POE();
    p = zeros(4,4,size(q,2));
    
    for i=1:size(q,2)
        p(:,:,i) = FKinSpace(M, Slist, q(:,i));
        p(:,:,i)
    end
     
    x = permute(p(1,end,:), [3 2 1]);
    y = permute(p(2,end,:), [3 2 1]);
    
    %close all
    figure(3)
    plot(x,y)
    hold on
    plot(x(1),y(1),'*')

end