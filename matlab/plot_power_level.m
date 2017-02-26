function plot_power_level( method_name, t, p, varargin )

global io

% get the legend name
switch method_name
    case 'brute_force_matlab'
        curr_leg = 'space-time-matlab';

    case 'solve_TD_diffusion'
        curr_leg = 'space-time';

    case 'solve_TD_diffusion_an_prec'
        curr_leg = 'space-time-ANALY';
    
    case 'solve_TD_diffusion_elim_prec'
        curr_leg = 'space-time-elim';
    
    case 'solve_IQS_diffusion_an_prec'
        curr_leg = 'IQS-an';
    
    case 'solve_IQS_diffusion_elim_prec'
        curr_leg = 'IQS-elim';
    
    case 'solve_IQS_diffusion_JFNK'
        curr_leg = 'IQS-JFNK';
    
    case 'solve_IQS_PC_diffusion_an_prec'
        curr_leg = 'IQS-PC-an';
    
    case 'solve_IQS_PC_diffusion_elim_prec'
        curr_leg = 'IQS-PC-elim';
    
    case 'solve_IQS_diffusion_td_prec'
        curr_leg = 'IQS-theta-prec';
    
    case 'solve_IQS_diffusion'
        curr_leg = 'IQS';
    
    case 'solve_PRKE'
        curr_leg = 'PRKE';
    
    case 'solve_PRKE_exact'
        curr_leg = 'PRKE exact';
    
    case 'solve_PRKE_QS'
        curr_leg = 'PRKE QS';
        
    otherwise
        whos method_name
        error('Method %s is UNKNOWN',method_name);
end

 
figure(io.figID); hold all;

% findobj(gcf,'Type','axes','Tag','legend')
a=get(legend(gca),'String');
if isempty(a)
    leg=char(curr_leg);
else
    leg=char(char(a),curr_leg);
end
plot(t,p,'o-');

% if IQS
nVarargs = length(varargin);
if nVarargs ==2
    t_iqs_fine = varargin{1};
    p_iqs_fine = varargin{2};
    plot(t_iqs_fine,p_iqs_fine,'c.-');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
end

legend(leg,'Location','Best');

hold off;


end
