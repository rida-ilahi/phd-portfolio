% Data Input
% File pairs: [Anterior, Posterior]
file_pairs = {
    'anterior1.csv', 'posterior1.csv';
    'anterior2.csv', 'posterior2.csv';
    'anterior3.csv', 'posterior3.csv';
};

base_path = ''; % Insert the base path to the folder where coordinates are stored
Fs = 8;           % Sampling rate (fps)
conv_factor = 6.1; % Pixels to micrometers
% linear
poly_order = 1;  % 1 = linear, 2 = quadratic, 3 = cubic

% Processes each cell
for f = 1:size(file_pairs, 1)
    
    % Load the posterior and anterior coordinates of the cell 
    ant_data = readtable(fullfile(base_path, file_pairs{f,1}));
    post_data = readtable(fullfile(base_path, file_pairs{f,2}));
    
    % The coordinates in the file are in px. To convert from px to um, use
    % the microscope conversion factor that depends on the objective of the
    % microscope. These are the coordinates in the lab frame.
    x_ant = ant_data.X / conv_factor;
    y_ant = ant_data.Y / conv_factor;
    x_post = post_data.X / conv_factor;
    y_post = post_data.Y / conv_factor;
    
    % Get the total time of the trajectory 
    N = length(x_ant);
    t = (0:N-1) / Fs;  % Time vector in seconds


    % Perform a linear fit on lab frame coordinates w.r.t time
    % y coordinates
    p_y_ant = polyfit(t, y_ant, poly_order);
    p_y_post = polyfit(t, y_post, poly_order);
    
    y_ant_fit = polyval(p_y_ant, t);
    y_post_fit = polyval(p_y_post, t);

    % x coordinates
    p_x_ant = polyfit(t, x_ant, poly_order);
    p_x_post = polyfit(t, x_post, poly_order);

    x_ant_fit = polyval(p_x_ant, t);
    x_post_fit = polyval(p_x_post, t);

    
    % Here, we change the coordinate system such that x-y coordinates align with the frame by frame movement of the cell instead of a global frame. 
    % This makes it possible to isolate the traverse movement from the forward motion, therefore enabling analysis of oscillatory parameters.  
    
    % Steps: 1) Calculate the midpoint of the x and y coordinates of the posterior and anterior ends of the
             %  cell.
           % 2) Express the lab frame coordinates relative to the
           %    midpoint of the coordinates
           % 3) Calculate the orientation of the lab frame coordinates (theta) with respect
           % to the lab frame coordinates.
           % 4) Use the orientation angle to define the rotation matrix.

    % Calculate the midpoint
    x_mid = (x_ant + x_post) / 2;
    y_mid = (y_ant + y_post) / 2;
    
    % Calculate cell axis orientation
    dx_cell = x_ant - x_post;
    dy_cell = y_ant - y_post;
    theta_cell = atan2(dy_cell, dx_cell); % This defines the orientation angle of cell axis w.r.t the lab frame!
    theta_smooth = movmean(unwrap(theta_cell), 5); % Prevents angle discontinutities (like jumping from -pi to pi)
                                                    % 5 is a smoothing
                                                    % average; it smooths
                                                    % the orientation angle over
                                                    % 5 frames.
                                                   
                                                
    % Calculates position of each coordinate relative to the midpoint
    x_ant_rel = x_ant - x_mid;
    y_ant_rel = y_ant - y_mid;
    x_post_rel = x_post - x_mid;
    y_post_rel = y_post - y_mid;
    
    % Initializing new body frame coordinates
    x_ant_body = zeros(N, 1);
    y_ant_body = zeros(N, 1);
    x_post_body = zeros(N, 1);
    y_post_body = zeros(N, 1);
    
    % Now, rotate points from the lab frame into the body frame
    % theta_smooth(i) is the body orientation in the lab frame, so to
    % express lab-frame coordinates in the body frame we apply the inverse
    % rotation. The inverse of a rotation by +theta is a rotation by
    % -theta.
    % In matrix form, this transformation would be: 
    % r_body = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)]*r_lab,
    % where [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)] is the inverted rotation matrix.

    for i = 1:N
        c = cos(-theta_smooth(i));
        s = sin(-theta_smooth(i));
        
        x_ant_body(i) = x_ant_rel(i) * c - y_ant_rel(i) * s;
        y_ant_body(i) = x_ant_rel(i) * s + y_ant_rel(i) * c;
        
        x_post_body(i) = x_post_rel(i) * c - y_post_rel(i) * s;
        y_post_body(i) = x_post_rel(i) * s + y_post_rel(i) * c;
    end
    
    % Fit a linear trend to the new body frame coordinates
    % polyfit finds the coefficients of the linear fit using a least-squares regression 
    % polyval(p, x) takes the coefficients p produced by polyfit and evaluates the polynomal at specified points x.

    p_y_ant_body = polyfit(t, y_ant_body, poly_order);
    p_y_post_body = polyfit(t, y_post_body, poly_order);
    
    y_ant_body_fit = polyval(p_y_ant_body, t);
    y_post_body_fit = polyval(p_y_post_body, t);

    p_x_ant_body = polyfit(t, x_ant_body, poly_order);
    p_x_post_body = polyfit(t, x_post_body, poly_order);

    x_ant_body_fit = polyval(p_x_ant_body, t);
    x_post_body_fit = polyval(p_x_post_body, t);
    
    % Subtracts linear trend to obtain oscillatory components
    y_ant_osc_body = y_ant_body - y_ant_body_fit;
    y_post_osc_body = y_post_body - y_post_body_fit;
    x_ant_osc_body = x_ant_body - x_ant_body_fit;
    x_post_osc_body = x_post_body - x_post_body_fit;

    %% Visualizing the coordinate system transformation with linear fits (3D plot with time in z)

    % Original coordinate system
    figure;
    plot3(x_post, y_post, t, 'b')
    hold on 
    plot3(x_post_fit, y_post_fit, t, 'b--')
    xlabel('x(\mum)');
    ylabel('y(\mum)'); % along the cell axis
    zlabel('t (s)'); % orthogonal to the cell axis
    box on

    % Transformed coordinate system
    figure;
    plot3(x_post_body, y_post_body, t, 'b')
    hold on 
    plot3(x_post_body_fit, y_post_body_fit, t, 'b--')
    xlabel('x(\mum)');
    ylabel('y(\mum)'); 
    zlabel('t (s)'); 
    box on

    %% This is the 2D projection of the posterior end oscillatory components of body frame coordinates 
    figure;
    yyaxis left
    plot(t, x_post_osc_body, 'b', 'LineWidth', 0.5)
    ylabel('x (\mum)')
    xlabel('t (s)')

    yyaxis right 
    plot(t, y_post_osc_body, 'r', 'LineWidth', 0.5)
    ylabel('y (\mum')

end