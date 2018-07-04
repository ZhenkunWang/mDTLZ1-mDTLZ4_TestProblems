function mop=mDTLZ(testname,pdim,odim)
%% The matlab code of mDTLZ problems.
%mDTLZ1-mDTLZ4 are four multi-objective optimization test problems
%proposed in the following paper:

% Wang Z, Ong Y S, and Ishibuchi H. "On Scalable Multiobjective Test
% Problems with Hardly Dominated Boundaries" IEEE Transactions on
% Evolutionary Computation, 2018, in press.

%% run for the test problems mdtlz serious.
%Get test multi-objective problems from a given name.The method is to get
%testing or benchmark problems for Multi-Objective Optimization.
%The test problem will be encapsulated in a structure,  which can be
%obtained by function get_structure('testmop').
The user gets the corresponding test problem, which is an instance of the class mop, by passing the problem name and optional dimension parameters.
mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
switch lower(testname)
    case 'mdtlz1'
        mop=mDTLZ1(mop,pdim,odim);
    case 'mdtlz2'
        mop=mDTLZ2(mop,pdim,odim);        
    case 'mdtlz3'
        mop=mDTLZ3(mop,pdim,odim);        
    case 'mdtlz4'
        mop=mDTLZ4(mop,pdim,odim);                                                     
    otherwise
        error('Undefined test problem name');
end
end
%%
%%%%%%%%%%FUNCTIONS%%%%%%
function p=mDTLZ1(p,pdim,odim)
   p.name='mDTLZ1';
   p.od=odim;
   p.pd=pdim;
   p.domain=[zeros(pdim,1) ones(pdim,1)];
   p.func=@evaluate;
   function y=evaluate(x)
      y=zeros(odim,1);
      g=zeros(odim,1);
      for j=1:odim
      xm=x(odim+j-1:odim:pdim,1);
      g(j)=100*(length(xm)+sum((xm-0.5).^2-cos(20*pi*(xm-0.5))));
      end
      y(1)=0.5*(1-prod(x(1:odim-1,1)))*(1+g(1));
      y(odim)=0.5*(x(1))*(1+g(odim));
      for i=2:odim-1
          y(i)=0.5*(1-prod(x(1:odim-i,1))*(1-x(odim-i+1)))*(1+g(i));
      end                        
   end
end
%%%%%%%
function p=mDTLZ2(p,pdim,odim)
   p.name='mDTLZ2';
   p.od=odim;
   p.pd=pdim;
   p.domain=[zeros(pdim,1) ones(pdim,1)];
   p.func=@evaluate;
   function y=evaluate(x)
      y=zeros(odim,1);
      g=zeros(odim,1);
      for j=1:odim
      xm=x(odim+j-1:odim:pdim,1);
      g(j)=sum((xm-0.5).^2);
      end
      y(1)=(1-prod(cos(pi*x(1:odim-1,1)./2)))*(1+g(1));
      y(odim)=(1-sin(pi*x(1)/2))*(1+g(odim));
      for i=2:odim-1
          y(i)=(1-prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2))*(1+g(i));
      end  
    end
end
%%%%%%%
function p=mDTLZ3(p,pdim,odim)
   p.name='mDTLZ3';
   p.od=odim;
   p.pd=pdim;
   p.domain=[zeros(pdim,1) ones(pdim,1)];
   p.func=@evaluate;
   function y=evaluate(x)
      y=zeros(odim,1);
      g=zeros(odim,1);
      for j=1:odim
      xm=x(odim+j-1:odim:pdim,1);
      g(j)=100*(length(xm)+sum((xm-0.5).^2-cos(20*pi*(xm-0.5))));
      end
      y(1)=(1-prod(cos(pi*x(1:odim-1,1)./2)))*(1+g(1));
      y(odim)=(1-sin(pi*x(1)/2))*(1+g(odim));
      for i=2:odim-1
          y(i)=(1-prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2))*(1+g(i));
      end       
    end
end
%%%%%%%
function p=mDTLZ4(p,pdim,odim)
   p.name='mDTLZ4';
   p.od=odim;
   p.pd=pdim;
   p.domain=[zeros(pdim,1) ones(pdim,1)];
   p.func=@evaluate;
   function y=evaluate(x)
      a=100;
      x=x.^a;
      y=zeros(odim,1);
      g=zeros(odim,1);
      for j=1:odim
      xm=x(odim+j-1:odim:pdim,1);
      g(j)=sum((xm-0.5).^2);
      end
      y(1)=(1-prod(cos(pi*x(1:odim-1,1)./2)))*(1+g(1));
      y(odim)=(1-sin(pi*x(1)/2))*(1+g(odim));
      for i=2:odim-1
          y(i)=(1-prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2))*(1+g(i));
      end  
   end
end
%%%%%%%

