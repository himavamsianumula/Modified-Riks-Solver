function project()
clc
clear all

% User inputs
th0=pi/3;
tol=1e-12;

%Newton's Method's inputs
incr_n=10000;
max_iter_n=60;
tol_n=1e-12;
dl=1.5e-2;

%Arc Length Method's inputs
incr_arc=2000;
max_iter_arc=10;
psi=1;
dll=6e-2;
tol_arc=1e-12;

[a_arc,al_arc]=riks(psi,dll,th0,incr_arc,max_iter_arc,tol_arc);
[a_t,al_t]=newton(dl,th0,incr_n,max_iter_n,tol_n);

%Plot for Newton's Method
figure(1)
plot(a_t,al_t,'x')
hold on
exx=linspace(0,2.6,501);
exy=((1./sqrt(b(exx,th0)) - 1.0).*(sin(th0)-exx));
plot(exx,exy)
title('Newton Method')
xlabel('Normalized Displacement')
xlim([0 2.5])
ylabel('Normalized force')
hold off

%Plot for Arc Length Method
figure(2)
plot(a_arc,al_arc,'o')
hold on
exx=linspace(0,2.6,501);
exy=((1./sqrt(b(exx,th0)) - 1.0).*(sin(th0)-exx));
plot(exx,exy)
title('Arc Length Method')
xlim([0 2.5])
xlabel('Normalized Displacement')
ylabel('Normalized force')
hold off

end

%% Newton's Method
function [a_t,al_t]=newton(dl,th0,incr,max_iter,tol)
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
                   disp([ 'Convergence cannot achieved within',max_iter,'iterations'])
                   disp('Program stops')
                   return
               end
           end
       end
    end
disp('The program completed successfully') 
end

%% Riks Method
function [a_t,al_t] = riks(psi,dll,th0,incr,max_iter,tol)
    
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
        if a(1)>2.5
            break;
        end
        da=zeros(n,1);
        dab=da;
        dat=da;
        dda1=da;
        dda2=da;
        dda=da;
        dalpha=da;
        dl=0.0;

        [df,dfinv]=dfcn(a+da,th0,al+dl);
        dat=dfinv'*iq;

        [ddl1,ddl2]=arc(psi,dll,da,dab,dat,dl,iq);

        dda1=dab+ddl1*dat;
        dda2=dab+ddl2*dat;

        %    det=det(df);
        det=df;

        if det*ddl1>0
            dda=dda1;
            ddl=ddl1;
        else
            dda=dda2;
            ddl=ddl2;
        end

        dalfa=da+dda;
        dlamda=dl+ddl;

        f=fcn((a+dalfa),th0,(al+dlamda));

        fcheck=sqrt(f'*f);

        if fcheck<tol
            a_t(count)=a;
            al_t(count)=al;
            a=a+dalfa;
            al=al+dlamda;
            
            count=count+1;
            a_t(count)=a;
            al_t(count)=al;
            
            dao=dalfa;
            dlo=dlamda;
            iloop=0;
        else
            iters=0;
            while fcheck>tol
                iters=iters+1;
                da=dalfa;
                dl=dlamda;

                f=fcn((a+da),th0,(al+dl));
                [df,dfinv]=dfcn((a+da),th0,(al+dl));

                dab=-dfinv'*f;
                dat=dfinv'*iq;

                [ddl1,ddl2]=arc(psi,dll,da,dab,dat,dl,iq);


                dda1=dab+ddl1*dat;
                dda2=dab+ddl2*dat;

                % 			det=det(df);
                det=df;

                daomag=(dao'*dao);

                if daomag==0.
                    if (dl+ddl1)*(det)>0
                        dda=dda1;
                        ddl=ddl1;
                    else
                        dda=dda2;
                        ddl=ddl2;
                    end

                else
                    aux1=(da+dda1)'*dao;
                    aux2=(da+dda2)'*dao;

                    aux3=dlamda*(dl+ddl1)*(iq'*iq);
                    aux4=dlamda*(dl+ddl2)*(iq'*iq);

                    dot1=aux1+(psi^2)*aux3;
                    dot2=aux2+(psi^2)*aux4;

                    if dot1>dot2
                        dda=dda1;
                        dd1=ddl1;
                    else
                        dda=dda2;
                        ddl=ddl2;
                    end

                end
                if ddl1==ddl2
                    dda=dda1;
                    ddl=ddl1;
                end

                dalfa=da+dda;
                dlamda=dl+ddl;

                f=fcn((a+dalfa),th0,(al+dlamda));
                fcheck=norm(f);

                if iters>max_iter
                    iters=max_iter+1;
                    break;
                end

            end

            if iters>max_iter
                disp(['Convergence cannot achieved within',max_iter,'iterations'])
                disp('Program stops')
                hold off
                return
            else
                a_t(count)=a;
                al_t(count)=al;
                a=a+dalfa;
                al=al+dlamda;
                
                count=count+1;
                a_t(count)=a;
                al_t(count)=al;
                dao=dalfa;
                dlo=dlamda;
            end
            plot(a,al,'x');
        end
    end
end
%% Helper Functions

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
	dfinv=inv(df);
end

function [ddl1,ddl2]=arc(psi,dll,da,dab,dat,dl,iq)
   
    c1=dat'*dat+(psi^2)*(iq'*iq);
    c2=2.0*(((da+dab)'*dat)+dl*(psi^2)*(iq'*iq));
    c3=(da+dab)'*(da+dab)+(dl^2)*(psi^2)*(iq'*iq)-(dll^2);
    
    if c2^2-4*c1*c3>0
        dls=roots([c1,c2,c3]);
        ddl1=dls(1);
        ddl2=dls(2);
    else
       ddl1=-c2/(2*c1);
       ddl2=-c2/(2*c1);
       disp('Poosible isssue in arc length')
    end
end
