if lambda_linear==2
    ind = find((Im==5|Im==6|Im==8|Im==9) & ym>mym+hydr_circ_thick); % find all markers from sediment 5, 6, 8, 9 where the vertical position below the surface line (mym) is greater than the y coordinate of the surface line (ym)
    lambdam(ind)= 0.4 + (ym(ind)-mym(ind)-hydr_circ_thick)./1000.*lambda_increase; % hydrostatic fluid P + (y-coordinate of surface elevation of ind - vertical position of markers ind) / (1000 m * 0.1 fluid pressure ratio), should increase fluid P by 0.1 per km
    ind= lambdam>lambda_bottom;    % where lambdam of ind is greater than lambda_bottom (0.9), make equal to 0.9
    lambdam(ind)=lambda_bottom;
end

if lambda_linear==3
    ind = find((Im==3|Im==4|Im==5|Im==6) & ym>mym+hydr_circ_thick);
    lambdam(ind)= 0.4 + (ym(ind)-mym(ind)-hydr_circ_thick)./(SediThick-hydr_circ_thick).*(lambda_bottom-0.4);% hydrostatic fluid P + (y-coordinate of surface elevation of ind - vertical position of markers ind) / (1000 m * 0.1 fluid pressure ratio), should increase fluid P by 0.1 per km
    ind = find((Im==3|Im==4|Im==5|Im==6) & ym>mym+SediThick);
    lambdam(ind)=lambda_bottom;      % lambda at sediment base
end

 ind = find(Im==13);
 lambdam(ind) = lambda_dry_OC;
 ind = find((Im==13) & Tm>200+273);
 lambdam(ind) = lambda_dry_OC + (lambda_hyd_OC-lambda_dry_OC).*(Tm(ind) - (200+273))./(100);
 ind = find((Im==13) & Tm>300+273);
 lambdam(ind)=lambda_hyd_OC;
%  ind = find((Im==12) & Tm>300+273);
%  lambdam(ind) = lambda_hyd_OC - (lambda_hyd_OC-lambda_dry_OC).*(Tm(ind) - (300+273))./(30);
%  ind = find((Im==12) & Tm>330+273);
%  lambdam(ind)=lambda_dry_OC;
