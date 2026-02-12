function [e] = fem_get_basis(p, quad_rule, etype);
%  fem_get_basis -- (FEM Tutorials)
%
%    This function generates a p^th order nodal Lagrange basis (default, p = 1), 
%    differentiates it diff times (default, diff=0), and returns the basis 
%    evaluated at all supplied xi-points.
%
%    Mandatory Inputs  --
%        p          - The nodal Lagrange polynomial order (integer > 0)
%
%        quad_rule - The quadrature rule order required
%
%        etype - Character variable 
%
%    Outputs --
%        e          - the basis object
%
%  by Dr. D Nordsletten
%  Jan 2015
%

  % set defaults ...

  % error checking ...
  fem_err(p <= 0, 'Polynomial order must be positive')
  fem_err(p - floor(p) ~= 0, 'Polynomial order must be an integer')

  fem_err(quad_rule <= 0, 'Quadrature order must be positive')
  fem_err(quad_rule - floor(quad_rule) ~= 0, 'Quadrature order must be an integer')

  fem_err(1-strcmpi(etype,'triangle')-strcmpi(etype,'quadrilateral'), 'Element type must be triangle or quadrilateral')
  
  e.etype = etype
  e.p = p;
  [e.gp, e.gw]=fem_get_quadrule_2D(quad_rule, etype); 
  [e.y]=fem_poly_2D(e.gp, etype, p);
  [e.dy(:,:,1)]=fem_poly_2D(e.gp, etype, p, [1 0]);
  [e.dy(:,:,2)]=fem_poly_2D(e.gp, etype, p, [0 1]);

  fem_check_basis(e)

