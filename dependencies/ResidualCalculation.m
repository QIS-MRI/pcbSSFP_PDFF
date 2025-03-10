clear res
TR = tr;
new_intense_mask = IntensityMask(profiles_lowres,mask_trs); % 

%%
profiles_lowres = profiles_lowres/abs(max(profiles_lowres(:)));
res = zeros(size(profiles_lowres,1), size(profiles_lowres,2), length(b0range));

for x = 1:size(profiles_lowres,1)
    for y = 1:size(profiles_lowres,2)

    p =  (squeeze(profiles_lowres(x, y, :)));
    
     k = 1;
    
        if new_intense_mask(x,y)>0
        
            p = ((squeeze(profiles_lowres(x, y, :))));
        
            if angle(mean(p))<0
                p = -1*p;
            end
        
            parfor ib0 = 1:length(b0range)
                b0 = b0range(ib0);
           
        
                [sol, quant,~,found_profile] = bSSFP_NonlinLS_Fit_WithPhase_Quant_MultComp(p, b0, b0+freqhz, amp, pc_step, TR, fa);
        
                t = norm(p-found_profile,2);
                res(x,y,ib0) = t;
            end
        else
        
        end
    end
        progress = x/size(profiles_lowres,1)*100
end


function concat_res = reimconcat(signal, dim)

    if nargin == 1
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = [signal_re;signal_im];
    else
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = cat(dim, signal_re, signal_im);
    end

end



function [Jm, F] = jf(x,J,f)

    Jm = J(x);
    F = f(x);

end

function e1 = getE1(a,b, alpha)
c = cosd(alpha);
e1 = (-b-c.*(b.*(a.^2))+a+(a.*c))./(a+(a.*c)-(b.*c)-(a.^2).*b);
end

function [sol, quant,res,found_profile] = bSSFP_NonlinLS_Fit_WithPhase_Quant_MultComp(p, phase, phase_fat, amps, pc_step, tr, fa)
    if nargin<3
    pc_step = 360/length(p);
    end
    pc = (0:pc_step:359)';
    pc = pc*pi/180;
    phase = phase*2*pi*(tr/1e3);
    pr = phase/2;

    phase_fat = phase_fat*2*pi*(tr/1e3);
    pr_fat = phase_fat*(1/2);

    ellipse_canon = @(x) (x(3)*(((1-x(1)*exp(1i*(-phase+pc))))./(1-x(2)*cos((-phase+pc)))))*exp(1i*pr);
    
    
    u0 = [ 0.99; 0.7; 1];
    u0_fat = [0.99; 0.7; 1];

    phase_kern = [cos(pr)*eye(length(pc)) -sin(pr)*eye(length(pc)); sin(pr)*eye(length(pc)) cos(pr)*eye(length(pc))];
    
   
    j1c = @(x) x(3)*(exp(1i*(-phase+pc))./(1-x(2)*cos(-phase+pc)));
    

    J1 = @(x) phase_kern*[real(j1c(x)); imag(j1c(x))];
    


    j2c = @(x) -x(3)* ( (1-x(1)*exp(1i*(-phase+pc))).*cos(-phase+pc) ./ ( ( 1-x(2)*cos(-phase+pc) ).^2 ) );
    
    J2 = @(x) phase_kern*[real(j2c(x)); imag(j2c(x))];
    

    j3c = @(x) -1*(((1-x(1)*exp(1i*(-phase+pc))))./(1-x(2)*cos(-phase+pc))); 
    

    J3 = @(x) phase_kern*[real(j3c(x)); imag(j3c(x))];
    
    for im = 1:length(phase_fat)

        phase_kern_fat{im} = [cos(pr_fat(im))*eye(length(pc)) -sin(pr_fat(im))*eye(length(pc)); ...
            sin(pr_fat(im))*eye(length(pc)) cos(pr_fat(im))*eye(length(pc))];
        j1c_fat = @(x) x(6)*(exp(1i*(-phase_fat(im)+pc))./(1-x(5)*cos(-phase_fat(im)+pc)));
        j2c_fat = @(x) -x(6)* ( (1-x(4)*exp(1i*(-phase_fat(im)+pc))).*cos(-phase_fat(im)+pc) ./ ...
            ( ( 1-x(5)*cos(-phase_fat(im)+pc) ).^2 ) );
        j3c_fat = @(x) -1*(((1-x(4)*exp(1i*(-phase_fat(im)+pc))))./(1-x(5)*cos(-phase_fat(im)+pc))); 
    
    
        J1_fat = @(x) phase_kern_fat{im}*[real(j1c_fat(x)); imag(j1c_fat(x))];
        J2_fat = @(x) phase_kern_fat{im}*[real(j2c_fat(x)); imag(j2c_fat(x))];
        J3_fat = @(x) phase_kern_fat{im}*[real(j3c_fat(x)); imag(j3c_fat(x))];

        ellipse_canon_fat = @(x) (x(6)*(((1-x(4)*exp(1i*(-phase_fat(im)+pc))))./(1-x(5)*cos((-phase_fat(im)+pc)))))*exp(1i*pr_fat(im));

        f = @(x) reimconcat( -ellipse_canon(x) - ellipse_canon_fat(x) ) + reimconcat(p);


        J = @(x) [J1(x)'; J2(x)'; J3(x)'; J1_fat(x)'; J2_fat(x)'; J3_fat(x)']';
        jf_inline{im} = @(x) jf(x,J,f);

    end



lambda = 1e-4;
u0 = [u0; u0_fat];
u = u0;
options = optimoptions("lsqlin");
options.Display = "off";
for i = 1:10
    Jm = zeros(length(reimconcat(p)), length(u));
    F =  zeros(length(reimconcat(p)), 1);
    for im = 1:length(phase_fat)
        jf_inline_p = jf_inline{im};
        [Jmd, Fd] = jf_inline_p(u);
        Jm = Jm+ amps(im)*Jmd;
        F = F + amps(im)*Fd;
    end
%     Jp = [Jm; sqrt(lambda)*eye(6)];
%     fp = [(Jm*u-F); sqrt(lambda)*u];
    b = Jm*u-F;
    A = Jm;
    try
    u = lsqlin(A,b,[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; ...
        -1 1 0 0 0 0; 0 0 0 -1 1 0; ...
        -1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],...
        [1 1 0.96 0.91 0 0 -0.82 -0.2 -0.82 -0.2], ...
        [],[],[],[],[],options);

    catch
        u = u0;
    end

end
    
    sol = u;
    found_profile = zeros(length(p),1);
    for im = 1:length(amps)
        ellipse_canon_fat = @(x) (x(6)*(((1-x(4)*exp(1i*(-phase_fat(im)+pc))))./(1-x(5)*cos((-phase_fat(im)+pc)))))*exp(1i*pr_fat(im));
        found_profile = found_profile + ( (1/length(phase_fat))*ellipse_canon(sol) + (amps(im))*ellipse_canon_fat(sol) );
    end
    res = norm(f(sol),2);
    a = abs(sol(1));
    q = abs(sol(2));
    M = abs(sol(3));
    a_fat = abs(sol(4));
    q_fat = abs(sol(5));
    M_fat = abs(sol(6));
    T2 = tr/log(1/a);
    E1 = abs(getE1(a,q,fa));
    T1 = abs(tr/log(1/E1));
    pd = (M/((1-E1)*sind(fa)/(1-E1*cosd(fa)-(sol(1)^2)*(E1-cosd(fa)))))...
        /(exp(-(tr/2)/T2));

    T2_fat = tr/log(1/a_fat);
    E1_fat = abs(getE1(a_fat,q_fat,fa));
    T1_fat = abs(tr/log(1/E1_fat));
    pd_fat = (M_fat/((1-E1_fat)*sind(fa)/(1-E1_fat*cosd(fa)-(sol(4)^2)*(E1_fat-cosd(fa)))))...
        /(exp(-(tr/2)/T2_fat));


    quant = [T1, T2, pd, T1_fat, T2_fat, pd_fat];
end