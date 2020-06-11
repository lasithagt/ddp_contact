%% Vector function for FK
function [pose_ret] = FK_SE3(Slist, M, q,i) 
    
    theta  = q;
%     pose_K = zeros(4,4,size(theta,2));
%     m = size(q,2);
    m = size(q,2);
    pose_K = zeros(4,4,m);
    

    
    for k=1:m
        pose_K(:,:,k)  = FKinSpace(M,Slist, theta(:,k));
        % pose_K(:,:,i)  = FK_kuka(theta(:,k));
    end
%     size(q,2)
%     i
    pose_ret = pose_K;
    
%     pose_ret = repmat(pose_K,1,1,size(q,2)/numel(i));
%     size(pose_ret);
    
   
end