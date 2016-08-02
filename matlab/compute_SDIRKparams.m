function [rk] = compute_SDIRKparams(time_integration)

switch time_integration
    case 1
        nstages=1;
        a=zeros(nstages,nstages);
        b=zeros(nstages,1); c=b;
        a(1,1)=1;
        c=sum(a,2);
        b(1)=1;
    case 2
        nstages=2;
        g=0.2928932188134524755991556378951510;
        a=[g, 0 ;...
          1-g, g];
        c=sum(a,2);
        b=a(nstages,:);
    case 3
        nstages=3;
        g=0.43586652150845899941601945119356;
        a=[g 0 0; ...
            ((1-g)/2) g 0;...
            (-(6*g^2-16*g+1)/4) ((6*g^2-20*g+5)/4) g];
        c=sum(a,2);
        b=a(nstages,:);
    otherwise
        error('time integration %d not yet implemented',time_integration);
end
rk.s=nstages;
rk.a=a; rk.b=b; rk.c=c;