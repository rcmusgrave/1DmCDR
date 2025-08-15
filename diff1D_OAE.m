function [ALK,DIC,F,TI,tiF]=diff1D_OAE(alk0,dic0,KAP,u,T,S,ti,z,params)
    %% Function for the integration of a 1D model for alkalinity and DIC
    %% INPUTS:
    %% initial conditions are alk0, dic0 [umol/kgSW]
    %% kap(nz,nt) is the diffusivity [m^2/s]
    %% u(nt) is the wind speed [m/s]
    %% T(nt) is the temperature at the surface [C]
    %% S(nt) is the salinity at the surface [psu]
    %% ti(nt) is the time in seconds for the integration
    %% z(nz) is the vertical coordinate [m]
    %% params is a structure with
    %%  params.dt is the integration timestep [s]
    %%  params.DT is the output timestep [s]
    %%  params.   contains all the CO2SYS parameters
    %% OUTPUTS:
    %% ALK, DIC are 2D arrays (z,TI) of alkalinity and DIC for integration period
    %% F is 1 1D array (ti) of flux at midpoints (i.e. ti(ii)-0.5*delt)
    %% RCM Apr 2024

    nt = length(ti);
    nz = length(z);
    dt = mean(diff(ti));
    dz = mean(diff(z));

    indt = floor(params.DT/dt); % output every indt iterations
    TI = ti(1:indt:end); % time vector for output fields
    NT = length(TI);

    ALK = NaN(nz,NT);
    DIC = NaN(nz,NT);
    F = NaN(1,nt-1);
    tiF = (ti(2:end)+ti(1:end-1))/2; % time grid for flux

    ALK(:,1) = alk0;
    DIC(:,1) = dic0;
    
    % % get pCO2 at surface to calculate flux
    % C_TADIC = CO2SYS(alk0(end),dic0(end),1,2,S(1),T(1),T(1),0,0,params.SI,params.PO4,params.NH4,params.H2S,params.pHSC,params.K1K2,...
    %     params.SO4,params.KF,params.BOR);
    % pco2 = C_TADIC(22);
    % [F_CO2, dpCO2]=FCO2(pco2,params.pCO2atm,T(1),S(1),u(1));
    % F(1) = F_CO2/(1000*24*3600); % [mol/m^2/s]

    %% parameters for Crank Nicolson 
    alp = dt/(2*dz^2); % 
    inc = 1;

    for ii = 2:nt

        rho = gsw_rho(S(ii),T(ii),0); % this is used to convert mol/kg -> mol/m^3
        kap = KAP(ii,:);
        
        kapHLF = NaN(1,nz+1); % kappa at half grid points
        kapHLF(2:end-1) = 0.5*(kap(2:end)+kap(1:end-1));
        kapHLF(1) = kap(1); % at i = 1/2
        kapHLF(end) = kap(end); % at i= N+1/2

        a = 1 - alp*(kapHLF(2:end)+kapHLF(1:end-1));
        RHS = spdiags([alp*kapHLF(2:end)' a' alp*kapHLF(1:end-1)'],[-1 0 1],nz,nz);     
        % note: The ordering of the diagonals is confusing
        RHS(1,1) = 1 - alp*kapHLF(2);
        RHS(end,end) = 1 - alp*kapHLF(end-1);
  
        b = 1 + alp*(kapHLF(2:end)+kapHLF(1:end-1));
        LHS = spdiags([-alp*kapHLF(2:end)' b' -alp*kapHLF(1:end-1)'],[-1 0 1],nz,nz);
        LHS(1,1) = 1 + alp*kapHLF(2);
        LHS(end,end) = 1 + alp*kapHLF(end-1);


        % THE BELOW IS WRONG

        % RHS = spdiags([alp*kapHLF(1:end-1)' a' alp*kapHLF(2:end)'],[-1 0 1],nz,nz); 
        % RHS(1,1) = 1 - alp*kapHLF(2);
        % RHS(end,end) = 1 - alp*kapHLF(end-1);
        % LHS = spdiags([-alp*kapHLF(1:end-1)' b' -alp*kapHLF(2:end)'],[-1 0 1],nz,nz);
        % LHS(1,1) = 1 + alp*(kapHLF(2));
        % LHS(end,end) = 1 + alp*kapHLF(end-1);

        LHSinv = inv(LHS);

        % get pCO2 at surface to calculate flux
        C_TADIC = CO2SYS(alk0(end),dic0(end),1,2,S(ii),T(ii),T(ii),0,0,params.SI,params.PO4,params.NH4,params.H2S,params.pHSC,params.K1K2,...
                        params.KSO4,params.KF,params.BOR);
        pco2 = C_TADIC(22);

        [F_CO2, dpCO2]=FCO2(pco2,params.pCO2atm,T(ii),S(ii),u(ii));
        F(ii-1) = F_CO2/(1000*24*3600); % [mol/m^2/s] this is at ti(ii)-0.5*delt

        dDICdz = F(ii-1)*1e6/kap(end)/rho;        % 

        RHS_alk = RHS*alk0;
        RHS_dic = RHS*dic0;
        RHS_dic(end) = RHS_dic(end) - 2*alp*kap(end)*dDICdz*dz; % umol/kgSW m/s
                                                                % RHS_dic(end) = RHS_dic(end) - 2*alp*F*dz; %% CHECK THIS

        alk1 = LHSinv*RHS_alk;               % ALK at next timestep
        dic1 = LHSinv*RHS_dic;               % DIC at next timestep

        alk0 = alk1;
        dic0 = dic1;

        if(mod(ii,indt)==0)          % save data every DT timesteps
            inc = inc+1;
            DIC(:,inc) = dic1;
            ALK(:,inc) = alk1;
        end
        
    end

end 