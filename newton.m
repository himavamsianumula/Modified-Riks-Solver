function [a_t,al_t,status]=newton(dl,th0,incr,max_iter,tol)
    n=1;
    iq=zeros(n,1);
    iq(1)=1;

    a=zeros(n,1);
    f=zeros(n,1);
    df=zeros(n,n);
    dfinv=zeros(n,n);

    dls=zeros(2,1);
    dao=zeros(n,1);
    al=0.0;
    
    a_t=0;
    al_t=0;
    
    count=1;
    
    for i=1:incr
       if a(1)>=2.5
           break;
       end
       da=zeros(n,1);
     
       a=a+da;
       al=al+dl;
       
       f=fcn((a),th0,(al));
       fcheck=sqrt(f'*f);
       
       if fcheck<tol
           iloop=0;
       else
           iters=0;
           while fcheck>tol
               iters=iters+1;
               [df,dfinv]=dfcn(a,th0,al);
               da=-1*(dfinv'*f);
               
               a=a+da;
               
               count=count+1;
               a_t(count)=a;
               al_t(count)=al;
               f=fcn(a,th0,al);
               fcheck=sqrt(f'*f);
               
               if iters>max_iter
%                    disp([ 'Convergence cannot achieved within',max_iter,'iterations'])
%                    disp('Program stops')
                        status=0;
                   return
               end
           end
       end
    end
 status=1;
% disp('The program completed successfully') 
end


function bb=b(x,y)
    bb=1.+x.^2.0-2.0.*x.*sin(y);
end

function f=fcn(x,y,z)
    bb=b(x(1),y);
	f(1)=(1./sqrt(bb)-1.0)*(sin(y)-x(1))-z;
end

function [df,dfinv]=dfcn(x,y,z)
    bb=b(x(1),y);
    df=zeros(1,1);
	df(1,1)=1-(1.-sin(y)^2.0)/(bb^1.5);
	dfinv=1/df;
end



