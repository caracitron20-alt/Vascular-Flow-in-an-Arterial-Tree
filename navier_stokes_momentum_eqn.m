function [Re] = navier_stokes_momentum_eqn(e,testsp, teste, vars);
% navier_stokes_eqn -- (FEM Tutorials)
%
% This function computes the element matrix for the 2D discrete Laplacian
%
% Inputs --
%
%        e       - the current element (integer)
%
%        testsp  - the space of test functions ...
%
%        teste   - the element / basis object for test functions 
%                  evaluated for the local element ...
%
%        vars    - a structure containing variables ...
%                  contains vars.[name]e - which contains basis functions 
%                                          evaluated for the local element
%
%
% Outputs --
%
%        Re - the element residual for Darcy-Navier-stokes equation ...
%
% by David Nordsletten
% Jan 2015
%

  % Problem parameters
 rho=1e-3;
 mu=4e-3;
  
  % Getting the local nodal values of velocity and pressure ...
  vl(:,1) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
  vl(:,2) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 2);
  pl(:,1) = vars.pres.u(vars.pres.t(e,:)); 

  % evaulating the weighted sum of velocity and pressure variables @ quadrature points ...
  v = [vars.vele.y(:,:) * vl(:,1) vars.vele.y(:,:) * vl(:,2)];
  p = vars.prese.y(:,:) * pl(:,1);


  % Velocity gradients
  dv1_dx1 = vars.vele.dy(:,:,1) * vl(:,1);
  dv1_dx2 = vars.vele.dy(:,:,2) * vl(:,1);
  dv2_dx1 = vars.vele.dy(:,:,1) * vl(:,2);
  dv2_dx2 = vars.vele.dy(:,:,2) * vl(:,2);

 
  % Getting local row / local column sizes ...
  ne=size(testsp.t,2);
  Re = zeros(testsp.dm * ne, 1);

  for i = 1:ne 
    % getting the local element residual indices ordered v_1^{n=1}, v_2^{n=1}, v_1^{n=2}, v_2^{n=2} ...
    vei = (testsp.dm * (i - 1) + 1):(testsp.dm * i);
    % Compute advection term: p(v·∇)v:
    advection_1 = rho * (v(:,1).* dv1_dx1+v(:,2).* dv1_dx2);
    advection_2 = rho * (v(:,1) .* dv2_dx1 + v(:,2) .* dv2_dx2);
    
    % computing ρ(v·∇v)·w + μ∇v:∇w - p∇·w
    Re(vei) = [dot(teste.gw, advection_1.* teste.y(:,i) ...
               + mu* (dv1_dx1.* teste.dy(:,i,1)+ dv1_dx2.* teste.dy(:,i,2)) ...
               - p.*teste.dy(:,i,1));
               dot(teste.gw, advection_2.* teste.y(:,i) ...
               + mu* (dv2_dx1.* teste.dy(:,i,1)+ dv2_dx2.* teste.dy(:,i,2)) ...
               - p.*teste.dy(:,i,2));];
  end


