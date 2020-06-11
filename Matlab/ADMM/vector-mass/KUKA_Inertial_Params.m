%% Global inertial parameters for the KUKA
function inertial_params = KUKA_Inertial_Params(lbr)

    frames = [{'iiwa_link_0'}, {'iiwa_link_1'}, {'iiwa_link_2'}, {'iiwa_link_3'}, {'iiwa_link_4'}, {'iiwa_link_5'},...
    {'iiwa_link_6'}, {'iiwa_link_7'}, {'iiwa_link_ee'}];
    for i = 2:7+1
        m(i-1) = lbr.Bodies{i}.Mass;
        I_ = lbr.Bodies{i}.Inertia;
        I_M = [I_(1) I_(6) I_(5);I_(6) I_(2) I_(4);I_(5) I_(4) I_(3)];
        mC(:,i-1) = lbr.Bodies{i}.CenterOfMass;
        I(:,:,i-1) = I_M;
    end
    inertial_params        = [m(1,:);mC(:,:);reshape([I(1,1,:);I(1,2,:);I(1,3,:);I(2,2,:);I(2,3,:);I(3,3,:)],[],7,1)];
    temp                   = [inertial_params(:);zeros(14,1)];
    inertial_params        = reshape(temp, 12, []);
    
    params    = inertial_params(:);
    params_fc = params(end-13:end);
    params    = params(1:70);
    params    = reshape(params, 10,[]);

    params    = [params;reshape(params_fc,2,[])];
    for i = 1 : 7
        T   = getTransform([0 0 0 0 0 0 0]',i,i+1);
        p   = T(1:3, end);  
        R   = T(1:3,1:3);
        params(2:4,i) = R'*(params(2:4,i) - p);

    end
     params(2:4,2) = params(2:4,2) + [0 0 0.2045]';
     params(2:4,4) = params(2:4,4) + [0 0 0.1845]';
     params(2:4,6) = params(2:4,6) + [0 0 0.081]';

%     for i = 1 : 7
%         T   = getTransform([0 0 0 0 0 0 0]',i,i+1);
%         p   = T(1:3, end);  
%         R   = T(1:3,1:3);
%         II     = [params(5,i) params(6,i) params(7,i);params(6,i) params(8,i) params(9,i);params(7,i) params(9,i) params(10,i)];
%         if (i==2 || i==4 || i==6)
%             temp_I = R*(II + 0*(params(1,i)*skew(params(2:4,i))'*skew(params(2:4,i))))*R';
%         else
%             temp_I = R*II*R';
%         end
%         params(5:10,i) = [temp_I(1,1) temp_I(1,2) temp_I(1,3) temp_I(2,2) temp_I(2,3) temp_I(3,3)];
%      end
%          
     % ---------------------------------------------------------------
    
    L_1xx = params(5,1); L_1xy = params(6,1); L_1xz = params(7,1); L_1yy = params(8,1); L_1yz = params(9,1); L_1zz = params(10,1); l_1x = params(2,1)*params(1,1); l_1y = params(3,1)*params(1,1); l_1z = params(4,1)*params(1,1); m_1 = params(1,1); fv_1 = params(11,1); fc_1 = params(12,1);
    L_2xx = params(5,2); L_2xy = params(6,2); L_2xz = params(7,2); L_2yy = params(8,2); L_2yz = params(9,2); L_2zz = params(10,2); l_2x = params(2,2)*params(1,2); l_2y = params(3,2)*params(1,2); l_2z = params(4,2)*params(1,2); m_2 = params(1,2); fv_2 = params(11,2); fc_2 = params(12,2);
    L_3xx = params(5,3); L_3xy = params(6,3); L_3xz = params(7,3); L_3yy = params(8,3); L_3yz = params(9,3); L_3zz = params(10,3); l_3x = params(2,3)*params(1,3); l_3y = params(3,3)*params(1,3); l_3z = params(4,3)*params(1,3); m_3 = params(1,3); fv_3 = params(11,3); fc_3 = params(12,3);
    L_4xx = params(5,4); L_4xy = params(6,4); L_4xz = params(7,4); L_4yy = params(8,4); L_4yz = params(9,4); L_4zz = params(10,4); l_4x = params(2,4)*params(1,4); l_4y = params(3,4)*params(1,4); l_4z = params(4,4)*params(1,4); m_4 = params(1,4); fv_4 = params(11,4); fc_4 = params(12,4);
    L_5xx = params(5,5); L_5xy = params(6,5); L_5xz = params(7,5); L_5yy = params(8,5); L_5yz = params(9,5); L_5zz = params(10,5); l_5x = params(2,5)*params(1,5); l_5y = params(3,5)*params(1,5); l_5z = params(4,5)*params(1,5); m_5 = params(1,5); fv_5 = params(11,5); fc_5 = params(12,5);
    L_6xx = params(5,6); L_6xy = params(6,6); L_6xz = params(7,6); L_6yy = params(8,6); L_6yz = params(9,6); L_6zz = params(10,6); l_6x = params(2,6)*params(1,6); l_6y = params(3,6)*params(1,6); l_6z = params(4,6)*params(1,6); m_6 = params(1,6); fv_6 = params(11,6); fc_6 = params(12,6);
    L_7xx = params(5,7); L_7xy = params(6,7); L_7xz = params(7,7); L_7yy = params(8,7); L_7yz = params(9,7); L_7zz = params(10,7); l_7x = params(2,7)*params(1,7); l_7y = params(3,7)*params(1,7); l_7z = params(4,7)*params(1,7); m_7 = params(1,7); fv_7 = params(11,7); fc_7 = params(12,7);

%     L_1xx = params(5,1); L_1xy = params(10,1); L_1xz = params(9,1); L_1yy = params(6,1); L_1yz = params(8,1); L_1zz = params(7,1); l_1x = params(2,1); l_1y = params(3,1); l_1z = params(4,1); m_1 = params(1,1); fv_1 = params(11,1); fc_1 = params(12,1);
% 
%     L_2xx = params(5,2); L_2xy = params(10,2); L_2xz = params(9,2); L_2yy = params(6,2); L_2yz = params(8,2); L_2zz = params(7,2); l_2x = params(2,2); l_2y = params(3,2); l_2z = params(4,2); m_2 = params(1,2); fv_2 = params(11,2); fc_2 = params(12,2);
% 
%     L_3xx = params(5,3); L_3xy = params(10,3); L_3xz = params(9,3); L_3yy = params(6,3); L_3yz = params(8,3); L_3zz = params(7,3); l_3x = params(2,3); l_3y = params(3,2); l_3z = params(4,3); m_3 = params(1,3); fv_3 = params(11,3); fc_3 = params(12,3);
% 
%     L_4xx = params(5,4); L_4xy = params(10,4); L_4xz = params(9,4); L_4yy = params(6,4); L_4yz = params(8,4); L_4zz = params(7,4); l_4x = params(2,4); l_4y = params(3,4); l_4z = params(4,4); m_4 = params(1,4); fv_4 = params(11,4); fc_4 = params(12,4);
% 
%     L_5xx = params(5,5); L_5xy = params(10,5); L_5xz = params(9,5); L_5yy = params(6,5); L_5yz = params(8,5); L_5zz = params(7,5); l_5x = params(2,5); l_5y = params(3,5); l_5z = params(4,5); m_5 = params(1,5); fv_5 = params(11,5); fc_5 = params(12,5);
% 
%     L_6xx = params(5,6); L_6xy = params(10,6); L_6xz = params(9,6); L_6yy = params(6,6); L_6yz = params(8,6); L_6zz = params(7,6); l_6x = params(2,6); l_6y = params(3,6); l_6z = params(4,6); m_6 = params(1,6); fv_6 = params(11,6); fc_6 = params(12,6);

%     L_7xx = params(5,7); L_7xy = params(10,7); L_7xz = params(9,7); L_7yy = params(6,7); L_7yz = params(8,7); L_7zz = params(7,7); l_7x = params(2,7); l_7y = params(3,7); l_7z = params(4,7); m_7 = params(1,7); fv_7 = params(11,7); fc_7 = params(12,7);

    inertial_params = [L_1xx L_1xy L_1xz L_1yy L_1yz L_1zz l_1x l_1y l_1z m_1 fv_1 fc_1...
        L_2xx L_2xy L_2xz L_2yy L_2yz L_2zz l_2x l_2y l_2z m_2 fv_2 fc_2...
        L_3xx L_3xy L_3xz L_3yy L_3yz L_3zz l_3x l_3y l_3z m_3 fv_3 fc_3...
        L_4xx L_4xy L_4xz L_4yy L_4yz L_4zz l_4x l_4y l_4z m_4 fv_4 fc_4...
        L_5xx L_5xy L_5xz L_5yy L_5yz L_5zz l_5x l_5y l_5z m_5 fv_5 fc_5...
        L_6xx L_6xy L_6xz L_6yy L_6yz L_6zz l_6x l_6y l_6z m_6 fv_6 fc_6...
        L_7xx L_7xy L_7xz L_7yy L_7yz L_7zz l_7x l_7y l_7z m_7 fv_7 fc_7];
    
    %% Outputs the transformation between given two frames in the first
    %  frame.
    function trans = getTransform(configuration,source, target)
        if nargin < 3
            source = target - 1; %if only target frame is specified, by default, output is with respect to the previous frame.
        end
%         if target > obj.NDOFs
%             warning('Target frame exceeds the number of degrees of the manipulator')
%         end
        
        source_ = frames{source+1};
        target_ = frames{target+1};
        trans  = lbr.getTransform(configuration,target_,source_);
        % R = trans(1:3,1:3); s = trans(1:3,end);
    end

end