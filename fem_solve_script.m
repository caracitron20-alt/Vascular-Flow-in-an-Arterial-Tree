%  AFEM solve 
%
%    This function assembles an FEM problem and solves it.  
%
%  by Cara Citron
%  Jan 2026
%
clear

  quad=10;
  etype = 'triangle';

  NoElems1D=10;
  h=[1 1] / NoElems1D;
  nonlinear_tol = 1e-6; % setting the nonlinear tolerance
  iter = 0; % initializing the interation counter
  residual_vector_norm = 1; % initializing the residual vector norm

  load('vesseltree3.mat')
% Making Domain function space ...
  Omega.e=fem_get_basis(Omega.p, quad, etype);
  Omega2.e=fem_get_basis(Omega2.p, quad, etype);
 

% Making variable function space ...
  vel=Omega2;
  vel.name='vel';
% Initializing the coefficients for the velocity field to zero (no nodes x number of unknowns per node)
  vel.u = zeros(vel.dm * size(vel.x,1), 1);
  n_vel = size(vel.u,1); 

% Making variable function space ...
  pres=Omega;
  pres.dm = 1;
  pres.name='pres';
% Initializing the coefficients for the pressure field to zero (no nodes x number of unknowns per node)
  pres.u = zeros(size(pres.x,1),1); 
  n_pres = size(pres.u,1);

while (iter == 0) || (residual_vector_norm > nonlinear_tol)
  % Making a FEM object storing my unknown variables (must do this every iteration) ...
    vars.vel = vel;
    vars.pres = pres;
    vars.Omega = Omega;

  % Making the residual vector and applying boundary conditions (velocity) ...
    R1 = fem_assemble_block_residual(@navier_stokes_momentum_eqn, Omega, vel, vars);
    for i = 1:size(vel.b,1)
      if(vel.b(i,5) == 2||vel.b(i,5) == 3||vel.b(i,5) == 4||vel.b(i,5) == 5) 
        continue
      end
        
     if(vel.b(i,5) == 1)
        for j=2:4
                nn = vel.dm * (vel.b(i,j) - 1);
                R1(nn+1:nn+2) = -[ -50 - vel.u(nn+1); 0 - vel.u(nn+2)];
        end
     end
      if(vel.b(i,5) == 6)
           for j=2:4
                nn = vel.dm * (vel.b(i,j) - 1);
                R1(nn+1:nn+2) = -[ 0 - vel.u(nn+1); 0 - vel.u(nn+2)];
           end
      end
    end
  % Making the residual vector and applying boundary conditions (on species B) ...
    R2 = fem_assemble_block_residual(@navier_stokes_mass_eqn, Omega, pres, vars);
    %R2(1) = 0;

  % Making the global residual vector + computing the norm ...
    R = [R1; R2];
    residual_vector_norm = norm(R,2);
    disp(['Current residual (iter=' num2str(iter) '): ' num2str(residual_vector_norm)])
    if(residual_vector_norm < nonlinear_tol)
      continue
    end

  % Creating the discrete matrix operators corresponding to operators a, b, and c
    disp(['        constructing the Jacobian blocks ...'])
    
    A  = fem_assemble_block_matrix_perturbation(@navier_stokes_momentum_eqn, Omega, vel, vel, vars); 
    B  = fem_assemble_block_matrix_perturbation(@navier_stokes_momentum_eqn, Omega, vel, pres, vars); 
    C  = fem_assemble_block_matrix_perturbation(@navier_stokes_mass_eqn, Omega, pres, vel, vars); 
    D = sparse(size(pres.u,1),size(pres.u,1));

  % Editing block matrices for conditions ...
    for i = 1:size(vel.b,1)
      if(vel.b(i,5) == 2||vel.b(i,5) == 3||vel.b(i,5) == 4||vel.b(i,5) == 5) 
        continue
      end
    
      if(vel.b(i,5) == 1 || vel.b(i,5) == 6)
        for j=2:4
                nn = vel.dm * (vel.b(i,j) - 1);
                A(nn+1:nn+2,:) = 0;
                A(nn+1:nn+2, nn+1:nn+2) = eye(2);
                B(nn+1:nn+2,:) = 0;
        end
     end

    end
    %C(1,:) = 0; D(1,1) = 1;

  % Composing the Jacobian from our Jacobian blocks ...
    disp(['        assembly of the Global Jacobian ...'])
    J = [ A B; C D];

  % Apply Newton Raphson ...
    disp(['        solving for the NR update ...'])
    U = J \ R;
    vel.u(1:n_vel) = vel.u(1:n_vel) - U(1:n_vel);
    pres.u(1:n_pres) = pres.u(1:n_pres) - U(n_vel+1:end);

  % Update the iteration counter ...
    iter = iter + 1;

    disp(['  ']) % skip a line so output is prettier
end

  % Plotting the velocity field ...
  quiver(vel.x(:,1), vel.x(:,2), vel.u(1:vel.dm:n_vel), vel.u(2:vel.dm:n_vel))

  % plotting the pressure field ...
  figure(2); trisurf(pres.t(:,1:3), pres.x(:,1), pres.x(:,2), pres.u,'Facecolor','interp','LineStyle','none')
  hold; trisurf(pres.t(:,1:3), pres.x(:,1), pres.x(:,2), pres.u,'Facecolor','interp','LineStyle','none')


%% Outflow rate trhouigh each terminal branch
% vel.x: node coordinates
% pres.u: nodal pressure 
% vel.u velocities 

%% Build element-to-edge map for outward normals
edge_to_elem = containers.Map('KeyType', 'char', 'ValueType', 'any');

for e = 1:size(vel.t, 1)
    elem_nodes = vel.t(e, 1:3);
    edges = [elem_nodes([1,2]); elem_nodes([2,3]); elem_nodes([3,1])];
    
    for i = 1:3
        edge_sorted = sort(edges(i,:));
        key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));
        
        if ~isKey(edge_to_elem, key)
            edge_to_elem(key) = [];
        end
        edge_to_elem(key) = [edge_to_elem(key); e];
    end
end

% display titles
disp('=========================================')
disp('Outflow  and Inflow rate through each terminal segment')

% Inlet flux for B1
inflow = 0;
inlet_edges = 0;

for i = 1:size(vel.b, 1)
    if vel.b(i, 5) ~= 1
        continue
    end
    
    n1 = vel.b(i, 2);
    n2 = vel.b(i, 3);
    
    x1 = vel.x(n1, :);
    x2 = vel.x(n2, :);
    L = norm(x2 - x1);
    
    u1 = vel.u(vel.dm*(n1-1)+1 : vel.dm*n1);
    u2 = vel.u(vel.dm*(n2-1)+1 : vel.dm*n2);
    u_avg = 0.5 * (u1 + u2);
    
    % Get outward normal
    edge = x2 - x1;
    normal1 = [edge(2), -edge(1)] / L;
    normal2 = [-edge(2), edge(1)] / L;
    
    edge_sorted = sort([n1, n2]);
    key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));
    
    if isKey(edge_to_elem, key)
        elem_list = edge_to_elem(key);
        if length(elem_list) == 1
            elem_idx = elem_list(1);
            elem_nodes = vel.t(elem_idx, 1:3);
            third_node = setdiff(elem_nodes, [n1, n2]);
            x_third = vel.x(third_node, :);
            edge_midpoint = 0.5 * (x1 + x2);
            to_third = x_third - edge_midpoint;
            
            if dot(normal1, to_third) < 0
                outward_normal = normal1;
            else
                outward_normal = normal2;
            end
        else
            outward_normal = normal1;
        end
    else
        outward_normal = normal1;
    end
    
    % For inlet, use inward normal (flip outward)
    inward_normal = -outward_normal;
    
    flux = dot(u_avg, inward_normal) * L;
    inflow = inflow + flux;
    inlet_edges = inlet_edges + 1;
end

disp(['Inlet (Γ_1): ' num2str(inflow, '%.6f') ' (' num2str(inlet_edges) ' edges)'])


% Outlet fluxes for BC 2,3,4 and 5 
boundary_labels = [2, 3, 4, 5];
outflow_rates = zeros(4, 1);

for k = 1:4
    bnd = boundary_labels(k);
    outflow = 0;
    edge_count = 0;
    
    edge_details = [];
    
    for i = 1:size(vel.b, 1)
        if vel.b(i, 5) ~= bnd
            continue
        end
        
        n1 = vel.b(i, 2);
        n2 = vel.b(i, 3);
        
        x1 = vel.x(n1, :);
        x2 = vel.x(n2, :);
        L = norm(x2 - x1);
        
        u1 = vel.u(vel.dm*(n1-1)+1 : vel.dm*n1);
        u2 = vel.u(vel.dm*(n2-1)+1 : vel.dm*n2);
        u_avg = 0.5 * (u1 + u2);
        
        % Get outward normal
        edge = x2 - x1;
        normal1 = [edge(2), -edge(1)] / L;
        normal2 = [-edge(2), edge(1)] / L;
        
        edge_sorted = sort([n1, n2]);
        key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));
        
        if isKey(edge_to_elem, key)
            elem_list = edge_to_elem(key);
            if length(elem_list) == 1
                elem_idx = elem_list(1);
                elem_nodes = vel.t(elem_idx, 1:3);
                third_node = setdiff(elem_nodes, [n1, n2]);
                x_third = vel.x(third_node, :);
                edge_midpoint = 0.5 * (x1 + x2);
                to_third = x_third - edge_midpoint;
                
                if dot(normal1, to_third) < 0
                    normal = normal1;
                else
                    normal = normal2;
                end
            else
                normal = normal1;
            end
        else
            normal = normal1;
        end
        
        flux = dot(u_avg, normal) * L;
        outflow = outflow + flux;
        edge_count = edge_count + 1;
        
        % Store for diagnostics
        edge_details = [edge_details; L, u_avg(1), u_avg(2), normal(1), normal(2), flux];
    end
    
    outflow_rates(k) = outflow;
    disp(['Outlet Γ_' num2str(bnd) ': ' num2str(outflow, '%.6f') ' (' num2str(edge_count) ' edges)'])
    
 end

total_outflow = sum(outflow_rates);

% Verify mass conservation and display
disp(' ')
disp('Mass Conservation:')
disp(['  Inflow:  ' num2str(inflow, '%.6f')])
disp(['  Outflow: ' num2str(total_outflow, '%.6f')])
disp(['  Error:   ' num2str(abs(inflow - total_outflow) / abs(inflow) * 100, '%.2f') '%'])

% Pressure calculation and display

disp('=========================================')
disp('Total pressure drop calculations')

% Initialise
inlet_pressure_sum = 0;
inlet_length = 0;

% Calculate Inlet Pressure 
for i = 1:size(vel.b, 1)
    if vel.b(i, 5) ~= 1
        continue
    end
    
    n1 = vel.b(i, 2);
    n2 = vel.b(i, 3);
    
    x1 = vel.x(n1, :);
    x2 = vel.x(n2, :);
    L = norm(x2 - x1);
    
    p1 = pres.u(n1);
    p2 = pres.u(n2);
    
    inlet_pressure_sum = inlet_pressure_sum + 0.5 * (p1 + p2) * L;
    inlet_length = inlet_length + L;
end

mean_inlet_pressure = inlet_pressure_sum / inlet_length;


% Calculate Outlet Pressure 

outlet_pressure_sum = 0;
outlet_length = 0;

for k = 1:4
    bnd = boundary_labels(k);
    
    for i = 1:size(vel.b, 1)
        if vel.b(i, 5) ~= bnd
            continue
        end
        
        n1 = vel.b(i, 2);
        n2 = vel.b(i, 3);
        
        x1 = vel.x(n1, :);
        x2 = vel.x(n2, :);
        L = norm(x2 - x1);
        
        p1 = pres.u(n1);
        p2 = pres.u(n2);
        
        outlet_pressure_sum = outlet_pressure_sum + 0.5 * (p1 + p2) * L;
        outlet_length = outlet_length + L;
    end
end

mean_outlet_pressure = outlet_pressure_sum / outlet_length;

disp(' ')
disp('Pressure:')
disp(['  Mean inlet pressure:  ' num2str(mean_inlet_pressure, '%.6f')])
disp(['  Mean outlet pressure: ' num2str(mean_outlet_pressure, '%.6f')])
disp(['  Total pressure drop:  ' num2str(mean_inlet_pressure - mean_outlet_pressure, '%.6f')])

%% Plotting for report
% Plot of convergence history (Figure 3)
figure(10);
clf;
iterations = 0:6;
semilogy(iterations, [141.42, 9.59, 3.88, 0.48, 8.8e-3, 4.7e-6, 4.1e-11], 'o-', 'LineWidth', 2, 'DisplayName', 'Mesh 0');
hold on;
semilogy(iterations(1:6), [200, 7.67, 2.28, 0.13, 5.3e-4, 1.6e-8], 's-', 'LineWidth', 2, 'DisplayName', 'Mesh 1');
semilogy(iterations(1:6), [278.39, 5.29, 1.61, 0.087, 2.7e-4, 3.8e-9], '^-', 'LineWidth', 2, 'DisplayName', 'Mesh 2');
semilogy(iterations(1:6), [400, 2.90, 0.99, 0.053, 1.7e-4, 3.3e-9], 'd-', 'LineWidth', 2, 'DisplayName', 'Mesh 3');
hold off;
xlabel('Iteration Number', 'FontSize', 12);
ylabel('Residual Norm', 'FontSize', 12);
title('Newton-Raphson Convergence History', 'FontSize', 14);
legend('Location', 'southwest');
grid off;
set(gca, 'FontSize', 11);

% Mesh Comparison (Mesh 0 vs Mesh 2) (Figure 4)
% Load mesh data for comparison
mesh_files = {'vesseltree0.mat', 'vesseltree2.mat'};
mesh_labels = {'Mesh 0', 'Mesh 2'};
mesh_data = cell(2, 1);

% Load and process both meshes
for m = 1:2
    mesh_temp = load(mesh_files{m});
    Omega_temp = mesh_temp.Omega;
    Omega2_temp = mesh_temp.Omega2;
    
    % Get basis functions
    Omega_temp.e = fem_get_basis(Omega_temp.p, quad, etype);
    Omega2_temp.e = fem_get_basis(Omega2_temp.p, quad, etype);
    
    % Setup velocity space (P2)
    vel_temp = Omega2_temp;
    vel_temp.name = 'vel';
    
    % Store mesh data
    mesh_data{m}.vel = vel_temp;
    mesh_data{m}.n_elements = size(vel_temp.t, 1);
    mesh_data{m}.n_dofs = vel_temp.dm * size(vel_temp.x, 1) + size(Omega_temp.x, 1);
end

% Create figure with two subplots
figure(14);
clf;
set(gcf, 'Color', 'w');

% Define boundary colors
colors = {'r', 'b', 'g', 'c', 'm', 'k'};
boundary_names = {'Γ_1 (Inlet)', 'Γ_2 (Outlet)', 'Γ_3 (Outlet)', 'Γ_4 (Outlet)', 'Γ_5 (Outlet)', 'Γ_6 (Wall)'};

% Plot both meshes
for m = 1:2
    subplot(1, 2, m);
    
    vel_temp = mesh_data{m}.vel;
    
    % Plot mesh triangulation
    triplot(vel_temp.t(:,1:3), vel_temp.x(:,1), vel_temp.x(:,2), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.3);
    hold on;
    
    % Highlight boundaries with colored lines
    for bnd = 1:6
        for i = 1:size(vel_temp.b, 1)
            if vel_temp.b(i, 5) == bnd
                n1 = vel_temp.b(i, 2);
                n2 = vel_temp.b(i, 3);
                plot([vel_temp.x(n1,1), vel_temp.x(n2,1)], ...
                     [vel_temp.x(n1,2), vel_temp.x(n2,2)], ...
                     colors{bnd}, 'LineWidth', 2.5);
            end
        end
    end
    
    title(sprintf('%s: %d elements, %d DOFs', ...
                  mesh_labels{m}, ...
                  mesh_data{m}.n_elements, ...
                  mesh_data{m}.n_dofs), ...
          'FontSize', 14);
    xlabel('x (mm)', 'FontSize', 12);
    ylabel('y (mm)', 'FontSize', 12);
    axis equal tight;
    grid on;
    hold off;
    set(gca, 'FontSize', 11);
end

% Add overall title
sgtitle('Mesh Refinement Comparison: Coarse vs Fine Resolution', 'FontSize', 16, 'FontWeight', 'bold');

% Add legend (create invisible plot for legend)
figure(15);
clf;
hold on;
legend_handles = [];
for bnd = 1:6
    legend_handles(bnd) = plot(NaN, NaN, colors{bnd}, 'LineWidth', 3);
end
legend(legend_handles, boundary_names, 'Location', 'best', 'FontSize', 12);
axis off;
title('Boundary Legend', 'FontSize', 14, 'FontWeight', 'bold');
hold off;

% Velocity Magnitude Contour (Figure 5)
figure(11);
clf;
vel_mag = sqrt(vel.u(1:vel.dm:end).^2 + vel.u(2:vel.dm:end).^2);
trisurf(vel.t(:,1:3), vel.x(:,1), vel.x(:,2), vel_mag, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
colorbar;
title('Velocity Magnitude Field (Mesh 3)', 'FontSize', 14);
xlabel('x (mm)', 'FontSize', 12);
ylabel('y (mm)', 'FontSize', 12);
axis equal tight;
set(gca, 'FontSize', 11);

% Pressure Field (Figure 6)
figure(12);
clf;
trisurf(pres.t(:,1:3), pres.x(:,1), pres.x(:,2), pres.u, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
colorbar;
title('Pressure Field Distribution (Mesh 3)', 'FontSize', 14);
xlabel('x (mm)', 'FontSize', 12);
ylabel('y (mm)', 'FontSize', 12);
axis equal tight;
set(gca, 'FontSize', 11);
