% pareto.m
% 
% The Matlab source codes to generate the PF and the PS of the test
% problems mDTLZ1-mDTLZ4 in the paper "Wang Z, Ong Y S, and Ishibuchi H.
% "On Scalable Multiobjective Test Problems with Hardly Dominated
% Boundaries" IEEE Transactions on Evolutionary Computation, 2018, in
% press."
% 
% Usage: [pf, ps] = pareto(problem_name, no_of_points, variable_dim)
% 
% Please refer to the report for more information.
% 
% History:
% v1 Jul.03 2018

function [pf, ps] = truePF_mDTLZ(name, no, dim)

    if nargin<3, dim = 3; end
    if nargin<2, no  = 500; end
    switch name
            case {'mDTLZ1'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = repmat(ps(2,:),[dim-2,1]).*repmat(ps(1,:),[dim-2,1]);             
            pf          = zeros(3,no);
            pf(1,:)     = 0.5*(1-ps(1,:).*ps(2,:));
            pf(2,:)     = 0.5*(1-ps(1,:).*(1-ps(2,:)));
            pf(3,:)     = 0.5*(ps(1,:));   
            clear s t;            
          case {'mDTLZ2','mDTLZ3','mDTLZ4'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = repmat(ps(2,:),[dim-2,1]).*repmat(ps(1,:),[dim-2,1]);             
            pf          = zeros(3,no);
            pf(1,:)     = 1-cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = 1-cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = 1-sin(0.5*pi*ps(1,:));   
            clear s t;            
    end
end
