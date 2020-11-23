% generate data points to learn a contact model
function contact_model_deep_learning

    % get samples
    
    
    % build a neural net.
    
    % train. try to use the GPU support.
    
    
    
    % test

end

function generate_samples(n_samples)
    
    %% generate random parameters
    m  = 1;
    mu = 0.4340;
    nu = 0.55;
    E  = 1e4;
    
    

    % generate data for the contact
    dyn_contact = contact_model_dyn(xd_e, xdd_e, f_ee, m, mu, nu, E);
    

end