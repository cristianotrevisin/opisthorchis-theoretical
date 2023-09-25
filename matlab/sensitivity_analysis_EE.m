function [UP,DOWN,fields] = sensitivity_analysis_EE(OCN,par,taboo)
    
    fields = fieldnames(par);
    UP = zeros(5,length(fields));
    DOWN = zeros(5,length(fields));

    % Run baseline
    [s,par] = build_setup(OCN,par,33*800*1000,'Unify',true);
    bsl = max(find_EE(par,par.c*s.H,s.H,s.S,s.KF,s.V));

    for pp = 1:length(fields)
        par.(fields{pp}) = par.(fields{pp})*1.2;
        [s,par] = build_setup(OCN,par,33*800*1000,'Unify',true);
        tmp = max(find_EE(par,par.c*s.H,s.H,s.S,s.KF,s.V));
        UP(:,pp) = (tmp-bsl)./bsl;
        
        par.(fields{pp}) = par.(fields{pp})/1.2*0.8;
        [s,par] = build_setup(OCN,par,33*800*1000,'Unify',true);
        tmp = max(find_EE(par,par.c*s.H,s.H,s.S,s.KF,s.V));
        DOWN(:,pp) = (tmp-bsl)./bsl;
        par.(fields{pp}) = par.(fields{pp})/0.8;
    end

    UP(abs(UP) < 1e-5)=0;
    DOWN(abs(DOWN)< 1e-5)=0;

    % Format labels for easier plotting
    fields = strrep(fields,'alpha','\alpha');
    fields = strrep(fields,'mu','\mu');
    fields = strrep(fields,'gamma','\gamma');
    fields = strrep(fields,'rho','\rho');
    fields = strrep(fields,'omega','\omega');
    fields = strrep(fields,'theta','\theta');
    fields = strrep(fields,'lambda','\lambda');
    fields = strrep(fields,'beta','\beta');
    fields = strrep(fields,'dF','d_F');
    fields = strrep(fields,'dS','d_S');


    for i = 1:length(taboo)
        %iDT = find([fields{:}]==[taboo(i)]);
        tmp = strfind(fields,taboo(i));
        iDT = find(not(cellfun('isempty',tmp)));
        UP(:,iDT) = [];
        DOWN(:,iDT) = [];
        fields(iDT) = [];
    end
    
    % Get order based on difference
    DIFF = abs(UP(1,:)-DOWN(1,:));
    [~,IDX] = sort(DIFF,'descend');
    UP_TEMP = UP;
    DOWN_TEMP = DOWN;
    fields_temp = fields;
    UP = zeros(5,length(fields_temp));
    DOWN = zeros(5,length(fields_temp));
    fields = cell(size(fields_temp));

    for rk = 1:length(DIFF)
        UP(:,rk) = UP_TEMP(:,IDX(rk));
        DOWN(:,rk) = DOWN_TEMP(:,IDX(rk));
        fields(rk) = fields_temp(IDX(rk));
    end



end