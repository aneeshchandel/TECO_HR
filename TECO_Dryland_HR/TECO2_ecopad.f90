! This version was built on TECO_spruce and modified by Enqing Hou for Sevilleta blue grama site
! By Enqing Hou in Nov. 2018

! *********************************************************
! Module notes
! Since the Sevilleta site do not have snow or ice. The processes of snow and ice were generally ignored.
! The current version does not work Methane or water table. Therefore, the two processes were also ignored.
! The current version added 

program TECO_MCMC
   
    implicit none  
!   Some initial settings       
    integer,parameter :: simu_n= 20000  ! number of simulation for DA
    integer,parameter :: rep_n= 100   ! number of forecasting
    
!   Initialize array setting  ************************************************************
!   Time series
    logical, parameter :: dosoilexp= .false.  ! for different soil types, including AMB
    logical, parameter :: doVcmaxexp= .false.
    logical, parameter :: dods0exp= .false.
    logical, parameter :: dogsc0exp= .false.
    integer year_length,yr_spinup
    integer, parameter :: day_length = 5114 ! 1096 ! determine size of blank array
    integer, parameter :: hour_length = day_length*24
    integer, parameter :: spinup_cycle = 0
    integer, parameter :: clim_var_n= 7  !  if(dosoilexp) clim_var_n= 11 else 7 ! number of climate input  ! for RCP scenarios
    
    integer, parameter :: year_length_ef = 18  ! year length for forecasting
    integer, parameter :: day_length_ef = 6574  ! day length for forecasting
    integer, parameter :: hour_length_ef = day_length_ef*24
    
!   number of parameters input from outside files    
    integer,parameter :: par_n=140    ! parameter in the main parameter file ("par.csv")
    integer,parameter :: par_stemp_n=20     ! 
    integer,parameter :: par_swater_n=6
    integer,parameter :: par_energy_n=13

    integer, parameter :: layern= 10  ! number of soil layers
    integer, parameter :: co2_n = 3   ! number of co2 scenarios
    
    integer,parameter :: swater_da_n=50     ! number of soil water parameters for data assimilation      
    integer, parameter :: cflux_obs_n=3   ! number of carbon fluxes
    integer, parameter :: swater_obs_n=6    ! numer of soil water observations

    integer, parameter :: simu_cflux_d_n =16     !   number of simulated carbon variables
    integer, parameter :: simu_swc_d_n = 14   ! number of simulated soil water variables
    

!   Switches *************************************************************************************    
    integer MCMC    !   0 for Simulation;  1 for Data assimilation;  2 for forecasting
    logical, parameter :: runlocal= .true.
    
    logical do_co2_da,do_swater_da
    logical, parameter :: do_soilt_da= .FALSE.
    logical, parameter :: do_snow_da= .False.
    logical, parameter :: do_methane_da= .False.
    logical, parameter :: do_soilphy= .true.
    logical, parameter :: do_snow   = .true.    ! if true, model simulate snow; else read snow from observation files

!   Switches for hourly write out
    logical, parameter :: wrt_cflux_d=.true.
    logical, parameter :: wrt_swater_d=.true.
    logical, parameter :: wrt_tsoil=.false.
    logical, parameter :: wrt_energy=.false.
    logical, parameter :: wrt_nflux=.false.
    logical, parameter :: wrt_swaterflux=.false.
    logical, parameter :: wrt_pheno=.false.
    logical, parameter :: wrt_cdist=.false.
!   Switches for yearly write out
    logical, parameter :: wrt_cflux_yr=.true.
    logical, parameter :: wrt_swater_yr=.true.
!   Switches for writing diagnose variables
    logical, parameter :: wrt_Aleaf=.false.
    logical, parameter :: wrt_acanop=.false.
    
    integer rcp
    real prescaler
    
!   Define directory --------------------------------------------------------------
    character(len=120) outdir,outdir_simu,outdir_forecast,outdir_da_cflux,outdir_da_swater    	
    character(len=60) outdir_da
!   1.1 Read TECO default parameters from pars.txt
    character(len=120) parafile
    real,dimension(par_n) :: par_main             
  
!   1.2 Read soil water parameters from pars_swater.txt, for soil_water subroutine
    character(len=120) parafile_swater
    real,dimension(par_swater_n) :: parval_swater
    real Crt,wpot50,condb

!   1.3 Read soil temperature parameters from pars_stemp.txt, for soil_temp subroutine
    character(len=120) parafile_stemp
    real,dimension(par_stemp_n) :: parval_stem    
    real shcap_snow,condu_snow,condu_b,depth_ex,albedo_snow,resht
    real fa,fsub,rho_snow,decay_m
    real shcap_soil,condu_soil,shcap_water,condu_water,shcap_air,condu_air
    real shcap_ice,condu_ice,latent_heat_fusion,ice_density
    
!   1.4 Read constant parameters from pars_constants.txt, for the whole model
    character(len=120) parafile_energy
    real,dimension(par_energy_n) :: parval_energy  
    real,dimension(3):: tauL,rhoL,rhoS
    real  emleaf,emsoil,gsc0,wleaf
    
!   1.5 Read initial soil conditions from soil_conditions.txt
    character(len=120) parafile_sconditions    
    real sthick(10),frlen(10),stemp(10),swc(10),waterv(10),ice(10)
    
!   1.7 Read initial site conditions from site_conditions.txt
    character(len=120) parafile_siteconditions
    real lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage
    
!   End of define parameters from files ********************************************
        
    
!   Blank parameters ***************************************************************
!   Forcing variables for TECO_simu

    character(len=150) climatefile,forcingdir
    real co2data(co2_n,hour_length)
    character(len=150) co2file
    integer first_year
    
!   Blank arrays for climate forcing used for simulation and data assimilation
    integer,dimension(hour_length):: year_seq,doy_seq,hour_seq
    real forcing_data(clim_var_n,hour_length)

    integer lines
 
!   Blank arrays for climate forcing used for ecological forecasting
    integer,dimension(hour_length_ef):: year_seq_ef,doy_seq_ef,hour_seq_ef
    real forcing_data_ef(clim_var_n,hour_length_ef)
    character(len=150) climatefile_ef   ! fixed name for TECO_simu
    integer lines_ef  
    
    
!   Blank arrays for observation data   
    real swater_obs(swater_obs_n,day_length),swater_obs_sd(swater_obs_n,day_length)  
    character(len=150) covfile,obsfile5,obsfile6       ! ..int 

    character(len=150) obs_nee_file
    integer obs_nee_d_date_temp(day_length),obs_nee_d_n
    integer,allocatable:: obs_nee_d_date(:)
    real obs_nee_d_temp(day_length)
    real,allocatable:: obs_nee_d(:)
    
!    aboveground carbon pool size
    character(len=150) obs_abc_file
    integer obs_abc_d_date_temp(day_length),obs_abc_d_n
    integer,allocatable:: obs_abc_d_date(:)
    real obs_abc_d_temp(day_length)
    real,allocatable:: obs_abc_d(:)
        
    character(len=150) obs_swc2p5_file
    integer obs_swc2p5_d_date_temp(day_length),obs_swc2p5_d_n
    integer,allocatable:: obs_swc2p5_d_date(:)
    real obs_swc2p5_d_temp(day_length)
    real,allocatable:: obs_swc2p5_d(:)
    
    character(len=150) obs_swc12p5_file
    integer obs_swc12p5_d_date_temp(day_length),obs_swc12p5_d_n
    integer,allocatable:: obs_swc12p5_d_date(:)
    real obs_swc12p5_d_temp(day_length)
    real,allocatable:: obs_swc12p5_d(:)
    
    character(len=150) obs_swc22p5_file
    integer obs_swc22p5_d_date_temp(day_length),obs_swc22p5_d_n
    integer,allocatable:: obs_swc22p5_d_date(:)
    real obs_swc22p5_d_temp(day_length)
    real,allocatable:: obs_swc22p5_d(:)
    
    character(len=150) obs_swc37p5_file
    integer obs_swc37p5_d_date_temp(day_length),obs_swc37p5_d_n
    integer,allocatable:: obs_swc37p5_d_date(:)
    real obs_swc37p5_d_temp(day_length)
    real,allocatable:: obs_swc37p5_d(:)
    
    character(len=150) obs_swc52p5_file
    integer obs_swc52p5_d_date_temp(day_length),obs_swc52p5_d_n
    integer,allocatable:: obs_swc52p5_d_date(:)
    real obs_swc52p5_d_temp(day_length)
    real,allocatable:: obs_swc52p5_d(:)
    
!   Blank arrays for simulation output    
    real simu_cflux_d(simu_cflux_d_n,day_length_ef),simu_swc_d(5,day_length)
    
!   Blank arrays for data assimilation   
    real parapost(200,10000)
    integer day_seq_spinup,seq,Pselect,wrtn
    real randv
    character(len=150) paraestfile
    integer,dimension(par_n):: da_check
    real,dimension(par_n) :: da_min,da_max
    character(len=par_n*15) indexstring
    
    
!   Others  **************************
    integer IDUM,upgraded,isimu     
    integer npara       ! Number of parameters to be estimated by DA
    real search_length,J_last,accR,J_temp
    real sdepth(10)
    
!   Blank arrays for regeneration of parameters __ Data assimilation
    real r

    real, allocatable :: coef_new(:), coef_old(:)
    real, allocatable :: coefmax(:),coefmin(:)
    integer,allocatable :: coefindex(:)
    integer k1,k2,rejet,paraflag,k3
    integer, parameter :: nc=100
    integer, parameter :: ncov=5000000
    character(len=150) outfile,MCMCargu,yrargu,dyargu
    character(len=150) Targu,CO2argu,comment
!   End of creat blank arrays ********************************************************
 
!   Blank arrays for initialization *********************************************************
    real swater_slw,swater_sld(10)
    real swater_slw_initial,swater_sld_initial(10)
    real Storage_initial
    real waterv_initial(10),ice_initial(10)
    integer rep,dylim
    character(len=150) my_fmt,fmt_cflux_da,fmt_swater_da
	
    character(len=150) yrlim_tem,dlim_tem
    character(len=21) cfluxfile
    character(len=22) swaterfile
    character(len=3) treatname,soiltype
    character(len=32) soilfiletemp
	
!   for simu_swc_d output    
    real wsc(10), water
    real covexist,randnum
    integer new,reject,tmp_up,i,j,k
    integer run,Jscaler

!    namelist
    character(len=50) nml_file
    namelist /control_pm/climatefile,forcingdir,parafile,parafile_siteconditions,parafile_sconditions, &
            parafile_stemp,parafile_energy, &
            outdir,outdir_simu,outdir_da,outdir_forecast, &
            first_year,year_length,yr_spinup, &
            MCMC,search_length,Jscaler,do_co2_da,do_swater_da, &
            obs_nee_file,obs_abc_file, &
            obs_swc2p5_file,obs_swc12p5_file, &
            obs_swc22p5_file, obs_swc37p5_file,obs_swc52p5_file
    
!    End of parameter defination  ***********************************************************************
    
!   Initialization ****************************************************************
!   Read in namelist
    print*,"Input name of namelist file:"
!    read(*,*) nml_file
    call getarg(1,nml_file)
!    nml_file = 'namelist_Ses2.nml'
    
    print*,nml_file
    open(1,file=trim(nml_file) )
    read(1,nml=control_pm)
    write(*,*)'climatefile=',climatefile
    write(*,*)'forcingdir=',forcingdir
    write(*,*)'parafile=',parafile
    write(*,*)'parafile_siteconditions=',parafile_siteconditions
    write(*,*)'parafile_sconditions=',parafile_sconditions
    write(*,*)'parafile_stemp=',parafile_stemp
    write(*,*)'parafile_energy=',parafile_energy
    write(*,*)'outdir=',outdir
    write(*,*)'outdir_simu=',outdir_simu
    write(*,*)'outdir_da=',outdir_da
    write(*,*)'outdir_forecast=',outdir_forecast
    write(*,*)'first_year=',first_year
    write(*,*)'year_length=',year_length
    write(*,*)'yr_spinup=',yr_spinup
    write(*,*)'MCMC=',MCMC
    write(*,*)'search_length=',search_length
    write(*,*)'Jscaler=',Jscaler
    write(*,*)'do_co2_da=',do_co2_da
    write(*,*)'do_swater_da=',do_swater_da
    write(*,*)'obs_nee_file=',obs_nee_file
    write(*,*)'obs_abc_file=',obs_abc_file
    write(*,*)'obs_swc2p5_file=',obs_swc2p5_file
    write(*,*)'obs_swc12p5_file=',obs_swc12p5_file
    write(*,*)'obs_swc22p5_file=',obs_swc22p5_file
    write(*,*)'obs_swc37p5_file=',obs_swc37p5_file
    write(*,*)'obs_swc52p5_file=',obs_swc52p5_file
    outdir = trim(outdir)
    
!   Initialize parameters for data assimilation
    J_last = 100000.00
    IDUM = 542
    upgraded=0
    new=0
    k3=0
    j=0
    parapost=0.0
    J_temp = 1
    
    swater_slw=1.0
    swater_sld=1.0

!   End of Initialization ****************************************************************
    
    call Getparameters_siteconditions(parafile_siteconditions, &
        &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage)

    
!   Define directory *******************************************************
           
	if(runlocal)then                
                outdir_da_cflux=trim(outdir_da)//'/cflux'
                outdir_da_swater=trim(outdir_da)//'/swater'   
                
                if(dosoilexp.or.doVcmaxexp.or.dogsc0exp.or.dods0exp)then
                    print*,"Input output directory:"
                    read(*,*) outdir_simu
                    print*,outdir_simu
                    print*,"Input climate file name:"
                    read(*,*) climatefile
                    print*,climatefile
                    print*,"Input precipitation scaler:"
                    read(*,*) prescaler
                    print*,prescaler
                    print*,"Input soil type name:"
                    read(*,*) soiltype
!                    print*,soiltype
                    parafile_sconditions=soilfiletemp //soiltype //".csv"
!                    print*,parafile_sconditions
                    print*,"Input RCP scenario:"
                    read(*,*) rcp
!                    print*,rcp
                    first_year=2018
                else
                    prescaler=1
!                    soiltype="AMB"
!                    parafile_sconditions=soilfiletemp // soiltype //".csv"
                    rcp=0
                endif
	 
	else		
		call getarg(1,parafile)
		call getarg(2,climatefile)
		call getarg(3,obs_nee_file)
		call getarg(4,outdir)
		call getarg(5,MCMCargu)
		read(MCMCargu,'(i1)') MCMC
		call getarg(6,parafile)
		call getarg(7,forcingdir)
		call getarg(8,yrlim_tem)
		call getarg(9,dlim_tem)
		call getarg(10,Targu)
		read(Targu,'(f9.3)') Ttreat
		call getarg(11,CO2argu) 
		read(CO2argu,'(f9.0)') CO2treat
		
		outdir_forecast="output"
		outdir_simu='output'
                outdir_da='output'
                outdir_da_cflux=outdir_da
                outdir_da_swater=outdir_da
		
	endif

    co2file="input/scenarios/co2conc.txt"

!   Switches for different scenarios ******************************************************

     
!   Read in inputs and parameters *************************************************************
    
!   1.1 Read parameters from file *************************************************************  
    call Getpara(parafile,par_n,    &
            &   par_main,da_check,da_min,da_max)
    npara=sum(da_check)
	print*, "Initial soil water contents are ", par_main(131:140)
    if(dods0exp)then
        print*,"Input Ds0:"
        read(*,*) par_main(4)
        print*,par_main(4)
    endif
            
    if(doVcmaxexp)then
        print*,"Input Vcmax,25:"
        read(*,*) par_main(5)
        print*,par_main(5)
    endif
              
!   1.2 Get parameters of soil water properties from pars_swater.txt, for soil_water subroutine 
!    parafile_swater='input/parameters/pars_swater.txt'
!    call Getparameters_swater(parafile_swater,crt,wlama,potcof,wpot50,condb,n_pot)  ! parameter related to soil water redistribution
!    parval_swater = (/crt,wlama,potcof,wpot50,condb,n_pot/)
   
!   1.3 Get parameters of soil thermal properties from pars_stemp.txt, for soil temperature subroutine
    
    call Getparameters_stemp(parafile_stemp,shcap_snow,condu_snow,condu_b,&  ! parameter related to snow
        &   depth_ex,albedo_snow,resht, &   ! parameter related to snow
        &   fa,fsub,rho_snow,decay_m,shcap_soil,condu_soil,         &   ! parameter related to snow
        &   shcap_water,condu_water,shcap_air,condu_air,&
        &   shcap_ice,condu_ice,latent_heat_fusion,ice_density)  
    parval_stem = (/shcap_snow,condu_snow,condu_b,& 
        &   depth_ex,albedo_snow,resht, &  
        &   fa,fsub,rho_snow,decay_m,shcap_soil,condu_soil,         &
        &   shcap_water,condu_water,shcap_air,condu_air,&
        &   shcap_ice,condu_ice,latent_heat_fusion,ice_density/)   

!    1.4 Get parameters of energy from pars_energy.txt
       
    call Getparameters_energy(parafile_energy,tauL,rhoL,rhos,emleaf,emsoil,wleaf,gsc0)
    parval_energy=(/tauL(1),rhoL(1),rhos(1),tauL(2),rhoL(2),rhos(2),tauL(3),rhoL(3),rhos(3),&
        &   emleaf,emsoil,wleaf,gsc0/)     
    
    if(dogsc0exp)then
        print*,"Input gsc0:"
        read(*,*) parval_energy(13)
        print*,parval_energy(13)
    endif    
        
!    1.5 Get initial soil conditions from  soil_conditions.txt

    call Get_soilconditions(parafile_sconditions,  &
    &   layern,sthick,frlen,stemp,waterv,ice)   
    swc=par_main(131:140)
    
!   End of reading parameters from files *****************************************
 
!   Read climate forcings *******************************************************
    if(rcp.gt.0)  call Getco2(co2file,co2_n,hour_length,co2data)
    
!   1 are climate data for simulation, 2 are climate data for forecasting   ..int
!   End of reading climate forcings  ***********************************************
     
    sdepth(1)=sthick(1)
    do i=2,layern
        sdepth(i)=sdepth(i-1)+sthick(i)
    enddo
    
!   Start main loop ********************************************************************************
    treatname=""
    if(dosoilexp)then
        treatname=soiltype
    endif
    
    if(doVcmaxexp.or.dogsc0exp.or.dods0exp) then
        print*,"Input file name:"
        read(*,*) treatname
        print*,"Treatment:",treatname
    endif

    call Getclimate(climatefile,hour_length,clim_var_n, &  ! input
            &   year_seq,doy_seq,hour_seq,forcing_data,lines) ! output
  

    if(MCMC.eq.0) then    ! Simulation
    print*,"Model simulation!!!"
    
    forcing_data(5,:)=forcing_data(5,:)*prescaler  
    
    
!   Write out simulated data  ---------------------------------
        if(wrt_cflux_yr) then
            if(rcp.eq.0)then
                cfluxfile="/simu_cflux_yr"//trim(treatname)//".txt"

                write(outfile,"(A120,A21)") trim(outdir_simu),cfluxfile
            else if (rcp.eq.85) then
                write(outfile,"(A120,A24)") trim(outdir_simu),"/simu_cflux_yr_rcp85.txt"
            else if  (rcp.eq.45)then
                write(outfile,"(A120,A24)") trim(outdir_simu),"/simu_cflux_yr_rcp45.txt"
            else if (rcp .eq.26)then
                write(outfile,"(A120,A24)") trim(outdir_simu),"/simu_cflux_yr_rcp26.txt"
            endif
                
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3301,file=outfile)
            write(3301,*) 'yr,GPP_yr,NPP_yr,NEE_yr,Resp_eco_yr,Resp_auto_yr, &
                    &   Resp_hetero_yr,NPP_L_yr,NPP_R_yr,LAI_yr'
        endif

        if(wrt_swater_yr) then
            if(rcp.eq.0)then
                swaterfile="/simu_swater_yr" //trim(treatname) // ".txt"
                write(outfile,"(A120,A22)") trim(outdir_simu),swaterfile
            else if (rcp.eq.85) then
                write(outfile,"(A120,A25)") trim(outdir_simu),"/simu_swater_yr_rcp85.txt"
            else if  (rcp.eq.45)then
                write(outfile,"(A120,A25)") trim(outdir_simu),"/simu_swater_yr_rcp45.txt"
            else if (rcp .eq.26)then
                write(outfile,"(A120,A25)") trim(outdir_simu),"/simu_swater_yr_rcp26.txt" 
            endif
            
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3302,file=outfile)
            write(3302,*) 'yr,transp_yr,evap_yr,runoff_yr,water_leach_yr,swater_slw_yr,swater_slw_gs, &
                    &   waterF1,waterF2,waterF3,waterF4,waterF5'
        endif
    
        if(wrt_cflux_d) then
            write(outfile,"(A120,A17)") trim(outdir_simu),"/simu_cflux_d.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3001,file=outfile)
            write(3001,*)'NEE,GPP,Reco,FoliageC'
        endif
        
        if(wrt_swater_d) then
            write(outfile,"(A120,A18)") trim(outdir_simu),"/simu_swater_d.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3002,file=outfile)
            write(3002,*) 'swc(1),swc(2),swc(3),swc(4),swc(5),swc(6),swc(7),swc(8),swc(9),swc(10)'
!            write(3002,*)'Unit,mm/h,mm/h,cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3,&
!                    &   cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3'
        endif
        
        if (wrt_tsoil) then
            write(outfile,"(A120,A15)") trim(outdir_simu),"/Simu_stemp.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3003,file=outfile)
            write(3003,*) 'day,Tsoil0,Tsoil1,Tsoil2,Tsoil3,Tsoil4,Tsoil5,Tsoil6,Tsoil7,&
                    & Tsoil8,Tsoil9,Tsoil10'
        endif
        
        if (wrt_energy) then
            write(outfile,"(A120,A16)") trim(outdir_simu),"/Simu_energy.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3004,file=outfile)
            write(3004,*) 'm,RsoilabT,Eabsorp_leaf,Esoil,Hsoil,G,Etransp,Hcan,Raero,gsc_l' 
        endif
        
        if (wrt_nflux) then
            write(outfile,"(A120,A15)") trim(outdir_simu),"/Simu_nflux.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3005,file=outfile)
            write(3005,*) 'day,Nuptake,Nminer,Nimmob,NresorpT,Nfix,Nleach,Nrunoff,Nvol,NdemandT,Ndeposit, &
                    &   NewplantN(1),NewplantN(2),NewplantN(3),MineralN,NSN,&
                    &   Npool(1),Npool(2),Npool(3),Npool(4),Npool(5),Npool(6),Npool(7),Npool(8), &
                    &   Nscal_vcmax,Nscal_autoresp'
        endif
        
        if (wrt_swaterflux) then
            write(outfile,"(A120,A21)") trim(outdir_simu),"/Simu_swater_flux.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3007,file=outfile)
            write(3007,*) 'm,water_free,runoff,water_leach,waterF(1),waterF(2),waterF(3),waterF(4),&
                             &  waterF(5),waterF(6),waterF(7),waterF(8),waterF(9),wpotent(1),wpotent(2),wpotent(3),&
                             &  wpotent(4),wpotent(5),wpotent(6),wpotent(7),wpotent(8),wpotent(9),wpotent(10)' 
!            write(3007,*) 'mm/h,mm/h,MPa,MPa,MPa,MPa,MPa,&
!                        & MPa,MPa,MPa,MPa,MPa,cm/h,cm/h,cm/h,cm/h,&
!                        & cm/h,cm/h,cm/h,cm/h,cm/h,&
!                        & water_fr,v/v,v/v,water_fr1,water_fr2,water_fr3,water_fr4,water_fr5,&
!                        & water_fr6,water_fr7,water_fr8,water_fr9,water_fr10' 
        endif
        
        if (wrt_pheno) then
            write(outfile,"(A120,A15)") trim(outdir_simu),"/Simu_pheno.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3100,file=outfile)
            write(3100,*) 'm,phenoset,onset,sene_count,&
                            &   GDD5,add,store,storage,accumulation'
        endif
        
        if (wrt_cdist) then
            write(outfile,"(A120,A24)") trim(outdir_simu),"/Simu_C_distribution.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3008,file=outfile)
            write(3008,*) 'GPP,NPPT,Resp_auto,add,store,Rmain,Rgrowth,Ccost_Nacq' 
        endif          
        
        if(wrt_Aleaf)then
            write(outfile,"(A120,A15)") trim(outdir_simu),"/Simu_Aleaf.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3200,file=outfile)
            write(3200,*) 'm,ng,swater_slw,Aleaf(1),Aleaf(2)' 
        endif
        
        if(wrt_acanop)then
            write(outfile,"(A120,A16)") trim(outdir_simu),"/Simu_Acanop.txt"
                outfile = trim(outfile)
                outfile = adjustl(outfile)
            open(3201,file=outfile)
            write(3201,*) 'm,swater_slw,Acanop,Acan(1),Acan(2),LAI' 
        endif
        
        write(*,*)'annual_summary:    year        LAI         gpp_yr         NPP_yr       ANPP_yr'
        call TECO_simu(MCMC,dosoilexp,do_soilphy,do_snow,rcp,clim_var_n,&  !   Switches
                &   year_length,day_length,hour_length,spinup_cycle,yr_spinup,& 
                &   year_seq,doy_seq,hour_seq,lines,day_seq_spinup,first_year,& !   time labels
                &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage, &   ! input: site conditions
                &   forcing_data,co2data,&  ! Forcing
                &   par_n,par_main,parval_swater,parval_stem,parval_energy,&  !   Grouped parameters
                &   sthick,frlen,layern,&   !    Site conditions
                &   swater_slw,swater_sld,swc,waterv,ice,wsc,& !   Soil water variables
                &   simu_cflux_d,   &   !   Simulation output
                &   simu_swc_d, & ! output
                &   wrt_cflux_d,wrt_swater_d,wrt_tsoil,wrt_energy,wrt_nflux,& ! switches: writing
                &   wrt_swaterflux,wrt_pheno,wrt_cdist,wrt_cflux_yr,wrt_swater_yr) ! switches: writing

        close(3001)
        close(3002)
        close(3003)
        close(3004)
         
        write(*,*)'run simulation'
        return   ! End of simulation (MCMC eq. 0)
        
    else if(MCMC.eq.1) then           ! Data assimilation 
        print*,"Data Assimilation!!!"
        
        swater_slw_initial = swater_slw
        swater_sld_initial = swater_sld
        Storage_initial = Storage

        write(outfile,"(A120,A12)") trim(outdir_da),"/Paraest.csv"
        outfile = trim(outfile)
        outfile = adjustl(outfile)
        open(4003,file=outfile)

        write(outfile,"(A120,A12)") trim(outdir_da),"/J value.csv"
        outfile = trim(outfile)
        outfile = adjustl(outfile)
        open(4004,file=outfile)

        print*,"Number of parameters for data assimilation:",npara
        allocate(coef_new(npara),coef_old(npara))
        allocate(coefindex(npara))
        allocate(coefmax(npara),coefmin(npara))

        do i=1,par_n
            if (da_check(i).eq. 1) then
                j=j+1
                coef_new(j)=par_main(i)        !define initial value for parameters, equal to the parameter file value
                coefindex(j)=i
                coefmin(j)=da_min(i)
                coefmax(j)=da_max(i)
            endif
        enddo
            
        write(4003,*)npara
        write(4003,*)(coefindex(i),i=1,npara)
        coef_old=coef_new
        
        if (do_co2_da)  then
            
!           Get observation data  
            if(obs_nee_file.ne.'')then
                call GetObs_d(obs_nee_file,day_length,obs_nee_d_n,obs_nee_d_date_temp,obs_nee_d_temp) 
                allocate(obs_nee_d_date(obs_nee_d_n))
                allocate(obs_nee_d(obs_nee_d_n))
                obs_nee_d_date=obs_nee_d_date_temp(1:obs_nee_d_n)
                obs_nee_d=obs_nee_d_temp(1:obs_nee_d_n)
            endif
            
            if(obs_abc_file.ne.'')then
                call GetObs_d(obs_abc_file,day_length,obs_abc_d_n,obs_abc_d_date_temp,obs_abc_d_temp) 
                allocate(obs_abc_d_date(obs_abc_d_n))
                allocate(obs_abc_d(obs_abc_d_n))
                obs_abc_d_date=obs_abc_d_date_temp(1:obs_abc_d_n)
                obs_abc_d=obs_abc_d_temp(1:obs_abc_d_n)
            endif
            
        endif

            
        if (do_swater_da) then

!           Get observation data
!           GetObs_d(obs_file,day_length,obs_d_n,obs_d_date_temp,obs_d_temp)
            if(obs_swc2p5_file.ne.'') then
                call GetObs_d(obs_swc2p5_file,day_length,obs_swc2p5_d_n,obs_swc2p5_d_date_temp,obs_swc2p5_d_temp)
                allocate(obs_swc2p5_d_date(obs_swc2p5_d_n))
                allocate(obs_swc2p5_d(obs_swc2p5_d_n))
                obs_swc2p5_d_date=obs_swc2p5_d_date_temp(1:obs_swc2p5_d_n)
                obs_swc2p5_d=obs_swc2p5_d_temp(1:obs_swc2p5_d_n)
            endif
            
            if(obs_swc12p5_file.ne.'') then
                call GetObs_d(obs_swc12p5_file,day_length,obs_swc12p5_d_n,obs_swc12p5_d_date_temp,obs_swc12p5_d_temp)
                allocate(obs_swc12p5_d_date(obs_swc12p5_d_n))
                allocate(obs_swc12p5_d(obs_swc12p5_d_n))
                obs_swc12p5_d_date=obs_swc12p5_d_date_temp(1:obs_swc12p5_d_n)
                obs_swc12p5_d=obs_swc12p5_d_temp(1:obs_swc12p5_d_n)
            endif
            
            if(obs_swc22p5_file.ne.'') then
                call GetObs_d(obs_swc22p5_file,day_length,obs_swc22p5_d_n,obs_swc22p5_d_date_temp,obs_swc22p5_d_temp)
                allocate(obs_swc22p5_d_date(obs_swc22p5_d_n))
                allocate(obs_swc22p5_d(obs_swc22p5_d_n))
                obs_swc22p5_d_date=obs_swc22p5_d_date_temp(1:obs_swc22p5_d_n)
                obs_swc22p5_d=obs_swc22p5_d_temp(1:obs_swc22p5_d_n)
            endif
            
            if(obs_swc37p5_file.ne.'') then
                call GetObs_d(obs_swc37p5_file,day_length,obs_swc37p5_d_n,obs_swc37p5_d_date_temp,obs_swc37p5_d_temp)
                allocate(obs_swc37p5_d_date(obs_swc37p5_d_n))
                allocate(obs_swc37p5_d(obs_swc37p5_d_n))
                obs_swc37p5_d_date=obs_swc37p5_d_date_temp(1:obs_swc37p5_d_n)
                obs_swc37p5_d=obs_swc37p5_d_temp(1:obs_swc37p5_d_n)
            endif
            
            if(obs_swc52p5_file.ne.'') then
                call GetObs_d(obs_swc52p5_file,day_length,obs_swc52p5_d_n,obs_swc52p5_d_date_temp,obs_swc52p5_d_temp)
                allocate(obs_swc52p5_d_date(obs_swc52p5_d_n))
                allocate(obs_swc52p5_d(obs_swc52p5_d_n))
                obs_swc52p5_d_date=obs_swc52p5_d_date_temp(1:obs_swc52p5_d_n)
                obs_swc52p5_d=obs_swc52p5_d_temp(1:obs_swc52p5_d_n)
            endif
            
        endif

        rejet = 0
        wrtn=0

    !    isimu=1
        do isimu=1,simu_n
    !       generate parameters
!            print*,"isimu:",isimu
            call gen_newcoef(coef_old,coefmax,coefmin,coef_new,search_length,npara)
            do k1=1,npara
                par_main(coefindex(k1))=coef_new(k1)
            enddo

            
            call TECO_simu(MCMC,dosoilexp,do_soilphy,do_snow,rcp,clim_var_n,&  !   Switches
                &   year_length,day_length,hour_length,spinup_cycle,yr_spinup, & 
                &   year_seq,doy_seq,hour_seq,lines,day_seq_spinup,first_year,& !   time labels
                &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage, &   ! input: site conditions
                &   forcing_data,co2data,&  ! Forcing
                &   par_n,par_main,parval_swater,parval_stem,parval_energy,&  !   Grouped parameters
                &   sthick,frlen,layern,&   !    Site conditions
                &   swater_slw,swater_sld,swc,waterv,ice,wsc,& !   Soil water variables
                &   simu_cflux_d,    &   !   Simulation output
                &   simu_swc_d, & ! output
                &   wrt_cflux_d,wrt_swater_d,wrt_tsoil,wrt_energy,wrt_nflux,& ! switches: writing
                &   wrt_swaterflux,wrt_pheno,wrt_cdist,wrt_cflux_yr,wrt_swater_yr) ! switches: writing

                
            tmp_up=upgraded

            call Costfunction(isimu,Jscaler,day_length,day_length_ef,do_co2_da,do_swater_da, &  ! input
                &   obs_nee_d_n,obs_nee_d_date,obs_nee_d,simu_cflux_d,   & ! input
                &   obs_abc_d_n,obs_abc_d_date,obs_abc_d,   & ! input
                &   obs_swc2p5_d_n,obs_swc2p5_d_date,obs_swc2p5_d,simu_swc_d,   & ! input
                &   obs_swc12p5_d_n,obs_swc12p5_d_date,obs_swc12p5_d,   & ! input
                &   obs_swc22p5_d_n,obs_swc22p5_d_date,obs_swc22p5_d,   & ! input
                &   obs_swc37p5_d_n,obs_swc37p5_d_date,obs_swc37p5_d,   & ! input
                &   obs_swc52p5_d_n,obs_swc52p5_d_date,obs_swc52p5_d,   & ! input
                &   J_temp,J_last,upgraded,accR) ! output
            
            if(upgraded.gt.tmp_up)then 
                new=new+1

                coef_old=coef_new
                call random_number(randv)
                write (fmt_cflux_da, '(a,i0,a,a)') '(F15.8,",",',npara-2,'(F15.8,",")','(F15.8))'
                write(4003,fmt_cflux_da) (coef_new(i),i=1,npara)  
                
                write(4004,40041) upgraded,J_last,accR
40041           format(i5,",",f10.4,",",f10.4)
                
!                if((upgraded.gt.10).and.(wrtn.lt.100).and.(randv.gt.0.0))then
                if((upgraded.gt.500).and.(wrtn.lt.100).and.(randv.lt.0.10)) then
                    wrtn=wrtn+1    
                    if(do_co2_da)then  

                        write(outfile,"(A120,A13,I3.3,A4)") trim(outdir_da_cflux), "/simu_cflux_d",wrtn,".txt"
                            outfile = trim(outfile)
                            outfile = adjustl(outfile)
                        open(4001,file=outfile)
                        write(4001,*)'NEE_d,GPP_d,Reco_d,FoliageC'
                        write(4001,*)'g/m2/d,g/m2/d,g/m2/d,g/m2'  
                        do i=1,day_length
                             write(4001,40011)(simu_cflux_d(j,i),j=1,4)
                        enddo           
    40011               format(3(f15.10,","),f15.10)
                    endif

                    if(do_swater_da)then

                        write(outfile,"(A120,A11,I3.3,A4)") trim(outdir_da_swater), "/simu_swc_d",wrtn,".txt"
                            outfile=trim(outfile)
                            outfile=adjustl(outfile)
                        open(4002,file=outfile) 
                        write(4002,*)'swc1,swc2,swc3,swc4,swc5' 
                        write(4002,*)'cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3'
                        do i=1,day_length
                            write(4002,40021) (simu_swc_d(j,i),j=1,5)
                        enddo           
    40021               format(4(f15.10,","),f15.10)       
                    endif

                endif  ! end of all true
                
            else
                reject=reject+1
            endif

            swater_slw = swater_slw_initial
            swater_sld = swater_sld_initial
            Storage=Storage_initial
            
        enddo ! End of loops from 1 to isimu

        close(4001)
        close(4002)
        close(4003)
        close(4004)

        deallocate(coef_new,coef_old)
        deallocate(coefindex)
        deallocate(coefmax,coefmin)
    !         
        write(*,*)'run MCMC'
   
    else                              ! (MCMC=2)   Forecasting
        print*,"Ecological forecasting!!!"
!    !   Update posterior parameters

        write(paraestfile,"(A120,A12)") trim(outdir_da),"/Paraest.csv"
        paraestfile = trim(paraestfile)
        paraestfile = adjustl(paraestfile)
        
        call Getparaest(paraestfile,par_n,parapost,seq,npara,indexstring)

        allocate(coefindex(npara))
        write (my_fmt, '(a,i0,a)') '(',npara,'I12)'
        read(indexstring,*) coefindex
       
        DO rep=1,rep_n

            CALL random_number(randnum)
            Pselect = int(seq/2+randnum*(seq-seq/2))
            do k1=1,npara
                par_main(coefindex(k1))=parapost(k1,Pselect)
            enddo
            
        !   Read generated climatic forcing
            !write(climatefile_ef,"(A120,A14,I3.3,A4)") trim(forcingdir),"/EMforcing_SEV",rep,".csv"
            write(climatefile_ef,"(A120)") trim(forcingdir)
            climatefile_ef=trim(climatefile_ef)
            climatefile_ef=adjustl(climatefile_ef)
!         print*,"test1:" ,climatefile_ef,hour_length_ef,clim_var_n
            call Getclimate(climatefile_ef,hour_length_ef,clim_var_n, &  ! input
                &   year_seq_ef,doy_seq_ef,hour_seq_ef,forcing_data_ef,lines_ef) ! output
   
            do k1=1,lines
                year_seq_ef(k1)=year_seq(k1)
                doy_seq_ef(k1)=doy_seq(k1)
                hour_seq_ef(k1)=hour_seq(k1)
                forcing_data_ef(1,k1)=forcing_data(1,k1)
                forcing_data_ef(2,k1)=forcing_data(2,k1)
                forcing_data_ef(3,k1)=forcing_data(3,k1)
                forcing_data_ef(4,k1)=forcing_data(4,k1)
                forcing_data_ef(5,k1)=forcing_data(5,k1)
                forcing_data_ef(6,k1)=forcing_data(6,k1)
                forcing_data_ef(7,k1)=forcing_data(7,k1)
            enddo

            call TECO_simu(MCMC,dosoilexp,do_soilphy,do_snow,rcp,clim_var_n, &  !   Switches
                    &   year_length_ef,day_length_ef,hour_length_ef,spinup_cycle,yr_spinup, & 
                    &   year_seq_ef,doy_seq_ef,hour_seq_ef,lines_ef,day_seq_spinup,first_year,& !   time labels 
                    &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage, &   ! input: site conditions
                    &   forcing_data_ef,co2data,&  ! Forcing
                    &   par_n,par_main,parval_swater,parval_stem,parval_energy,&  !   Grouped parameters
                    &   sthick,frlen,layern,&   !    Site conditions
                    &   swater_slw,swater_sld,swc,waterv,ice,wsc, & !   Soil water variables
                    &   simu_cflux_d,   &   !   Simulation output
                    &   simu_swc_d, & ! output
                    &   wrt_cflux_d,wrt_swater_d,wrt_tsoil,wrt_energy,wrt_nflux,& ! switches: writing
                    &   wrt_swaterflux,wrt_pheno,wrt_cdist,wrt_cflux_yr,wrt_swater_yr) ! switches: writing
                         
            write(*,*)'run forecasting',rep

            write(outfile,"(A120,A13,I3.3,A4)") trim(outdir_forecast), "/simu_cflux_d",rep,".txt"
            outfile=trim(outfile)
            outfile=adjustl(outfile)
            open(5001,file=outfile)
            write(5001,*)'NEE_d,GPP_d,Resp_eco_d,FoliageC'
            write(5001,*)'g/m2/d,g/m2/d,g/m2/d,g/m2'  
            do i=1,day_length_ef
                write(5001,50011) i, (simu_cflux_d(j,i),j=1,4)
            enddo           
50011       format((I5,","), 3(f15.10,","),f15.10)
            close(5001)
            
            write(outfile,"(A120,A11,I3.3,A4)") trim(outdir_forecast), "/simu_swc_d",rep,".txt"
            outfile=trim(outfile)
            outfile=adjustl(outfile)
            open(5002,file=outfile) 
            write(5002,*)'swc1,swc2,swc3,swc4,swc5' 
            write(5002,*)'cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3,cm3/cm3'
            do i=1,day_length
                write(5002,50021) (simu_swc_d(j,i),j=1,5)
            enddo           
50021       format(4(f15.10,","),f15.10)    
            close(5002)

        enddo ! END of replications of forecasting
        
    endif ! End of if(MCMC)

end  !    End of TECO_MCMC
!   ====================================================================================================


!  =====================================================================================================
subroutine TECO_simu(MCMC,dosoilexp,do_soilphy,do_snow,rcp,clim_var_n,&  !  input: Switches
    &   year_length,day_length,hour_length,spinup_cycle,yr_spinup, & ! input: timeline
    &   year_seq,doy_seq,hour_seq,lines,day_seq_spinup,first_year,& !  input:   time labels 
    &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage, &   ! input: site conditions
    &   forcing_data,co2data,&  !  input: Forcing
    &   par_n,par_main,parval_swater,parval_stem,parval_energy,&  !  input:   Grouped parameters
    &   sthick,frlen,layern,&   !   input:   Site conditions
    &   swater_slw,swater_sld,swc,waterv,ice,wsc, & !  input:    Soil water variables
    &   simu_cflux_d,    &   !  output:   Simulation output
    &   simu_swc_d, & ! output
    &   wrt_cflux_d,wrt_swater_d,wrt_tsoil,wrt_energy,wrt_nflux,& ! switches: writing
    &   wrt_swaterflux,wrt_pheno,wrt_cdist,wrt_cflux_yr,wrt_swater_yr) ! switches: writing

    implicit none
    
!   Switches
    logical,parameter::tec_diagn=.false.
	logical,parameter::rechargeDS=.false.
    logical do_snow,do_soilphy,dosoilexp
    logical wrt_cflux_d,wrt_swater_d,wrt_tsoil,wrt_energy,wrt_swaterflux
    logical wrt_cflux_yr,wrt_swater_yr
    logical wrt_pheno,wrt_nflux,wrt_cdist
    integer,parameter:: vegetype=1
    integer rcp
    
    
!   Settings
    integer layern
    
!   1.1 Read TECO default parameters from pars.txt
    integer par_n
    real,dimension(par_n) :: par_main 
    real lat,longi,wsmax(10),wsmin(10)                 
    real rdepth,SLA,stom_n,Ds0,Vcmax
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,gcostpro,mresp20(3),Q10paccrate(3),basew4sresp
    
!   1.2 Read soil water parameters from pars_swater.txt, for soil_water subroutine     
    integer,parameter :: par_swater_n=6
    real,dimension(par_swater_n) :: parval_swater
    real watersat(10),waterres(10),condsat(10),wlama(10),potcof(10),n_pot(10)

!   1.3 Read soil temperature parameters from pars_stemp.txt, for soil_temp subroutine    
    integer,parameter :: par_stemp_n=20
    real,dimension(par_stemp_n) :: parval_stem

!   1.4 Read parameters of energy from pars_energy.txt, for the whole model    
    real parval_energy(13)
    real,dimension(3):: tauL,rhoL,rhoS
    real  pi,emleaf,emsoil,Rconst,sigma,cpair,Patm,Temp_ref,H2OLv0,AirMa,H2OMw,chi,Dheat
    real  gsc0,wleaf,Kc_ref,conKo0,Ekc,Eko,o2ci
    real  gam0,gam1,gam2,Jmax_n,Vcmax_n,stress_scal,gdd_scal,LAI_scal
    
    
    
!   1.# Climate forcing
    integer clim_var_n            ! 9 for Duke forest FACE
    integer, parameter :: co2_n=3            ! 9 for Duke forest FACE
    real, parameter:: times_storage_use=720.   ! Initial: 720 hours, 30 days
    integer  lines,day_of_year,MCMC
    integer,dimension(hour_length):: year_seq,doy_seq,hour_seq
    real forcing_data(clim_var_n,hour_length)
    real co2ca,wind,esat,Dair,rain,radsol,Tair,RH,tair_dmean
    real TairK
    
    real co2data(co2_n,hour_length)

!   1.7 Read initial site conditions from site_conditions.txt
    real Ttreat,CO2treat,Ndeposit,Nfert,Storage
    
!   1.* Write out options
    real simu_cflux_d(16,day_length),simu_swc_d(5,day_length)
   
  
!   Phenology ----------------------------------------------------
    integer pheno,phenoset,onset !flag of phenological stage
    
!   C pool and fluxes -------------------------------------------
    real TauC(8),Cpool(8),OutC(8),Rh_pools(5)
    real NSNcrt1,NSNmin,add             ! none structural carbon pool
    real gpp,NPPT,Acanop
    real GLmx,Gsmx,GRmx
    real rootsap,stemsap
    
    real,dimension(5):: Gaussx,Gaussw,Gaussw_cum 
    real LAI
    real litterfall,rootfall
    real GDD5,accumulation,stor_use,store
    integer sene_count
    real RaL,RaS,RaR  !allocation to respiration
    real Ccost_Nacq   ! respirations
    real Resp_auto,Resp_hetero,Resp_eco,NEE !respirations
    real RmLeaf,RmStem,RmRoot          ! maintanence respiration
    real RgLeaf,RgStem,RgRoot          ! growth respiration

    real NPP(3),Growth(3)
    real gsc_l
    
!   Plant respiration ------------------------------------------
    
    
!   Soil_water subroutine --------------------------------------
    real swater_slw,swater_sld(10),relsat(10),wsc(10),swc(10),waterv(10),ice(10)
    real wsmin_fr,swc_fr,waterF(10),wpotent(10)
    real runoff,water_free,water_leach 
    real Rsoil
    
!   Soil temperature properties --------------------------------
    real sftmp,raero
    real stemp(10)
    real Tsoil  
    real Tleaf,Troot
    real Esoil_res
  
!   Model settings -----------------------------------------------
    real sthick(10),FRLEN(10)   ! wsc is the output from soil water module

!   Water fluxes ---------------------------------------------- 
    real evap,transp
    
!   Energy fluxes ---------------------------------------------
    real G,Esoil,Hsoil,Eabsorp_leaf,Etransp,Hcan
    real Radabv(2)
    real Qcan(3,2)
    
!   Snow parameters ----------------------------------------------
    real thd_snow_depth,water_table_depth,rain_d,zwt,Tsnow   
    real snow_depth,dcount
    real ice_tw
    real decay_m,fa,fsub,melt,rho_snow
    
!   Canopy subroutine -------------------------------------------
    real fbeam
    real RsoilabT,Rsoilab(3)
    real rhocp,Cmolar,H2OLv,psyc,slope
    
!   plant_growth subroutine --------------------------------------
    real growthT
    
!   plant_respiration subroutine --------------------------------
    real RmainT,Rmain(3),RgrowthT,Rgrowth(3)
    
!   For Eco_N subroutine ----------------------------------------
    real CNini(8),CN(8)
    real NSN,MineralN,Npool(8),OutN(8)
    real NdemandT,Ndemand,NresorpT,NimmobT
    real Ndeficit
    real Nuptake,Nfix,Nminer,Nleach,Nrunoff,Nvol,Nimmob,Nresorp,NewplantN(3)
    real Nscal_vcmax,Nscal_autoresp
    
!   Numbers, date, time ------------------------------------------------------
    integer doy,hour
    integer year,yr,days,i,j,m,yr_spinup,iyr,idays,day_n
    integer day_seq_spinup,day_onecycle,day_finalcycle,spinup_cycle
    integer hour_of_day,year_length,day_length,hour_length
    integer first_year
    integer start_growing_season,end_growting_season
    
!   Yearly output ------------------------------------------------
    real GPP_yr,NPP_yr,ANPP_yr,NEE_yr
    real Resp_eco_yr,Resp_auto_yr,Resp_hetero_yr,NPP_L_yr,NPP_R_yr,LAI_yr
    
    real transp_yr,evap_yr,runoff_yr,water_leach_yr,swater_slw_yr,swater_slw_gs
    real waterF_yr(10)
    
!   Daily output -------------------------------------------------
    real NEE_d,GPP_d,Resp_eco_d,foliageC,swc_d(5)

    
!   End of parameter defination -------------------------------------------------
    
    
    
!   Constants ----------------------------------------------------
    Rconst=8.314 ! Universal gas constant; unit of J/mol
    Patm=101325     ! Atmospheric pressure; unit of Pa
    cpair=1010      ! Heat capacity of air/air specific heat capacity; unit of J/kg/K
    H2OMw=1.80E-02  ! Molar mass of H2O; unit: kg/mol
    airMa=2.90E-02  ! Molar mass of air; unit: kg/mol
    sigma=5.67E-08  ! Steffan Boltzman constant; unit: W/m2/K4
    H2OLv0=2.50E+06 ! Latent heat of H2O; unit: J/kg

    
!   *******************************************88888888888888888888888
    
!    if (MCMC.eq.0) then
!        print*,"Simulation start"
!        print*,"yr_spinup:",yr_spinup,"          year_length:", year_length
!    endif
    
    
    rdepth=par_main(1)
    SLA=par_main(2)
    stom_n=par_main(3)
    Ds0=par_main(4)
    Vcmax=par_main(5)
    tauC=par_main(6:13)*8760
    gddonset=par_main(14)
    Q10=par_main(15)
    gcostpro=par_main(16)
    mresp20=par_main(17:19)
    Q10paccrate=par_main(20:22)
!    par_main=par_main(26:33)
    basew4sresp=par_main(34)
    wsmax=par_main(35:44)
    wsmin=par_main(45:54)
    watersat=par_main(55:64)
    waterres=par_main(65:74)
    condsat=par_main(75:84)
    wlama=par_main(85:94)
    potcof=par_main(95:104)
    n_pot=par_main(105:114)
    Cpool=par_main(115:122)
    CNini=par_main(123:130)
    
    decay_m=parval_stem(15)
    fa=parval_stem(12)
    fsub=parval_stem(13)
    rho_snow=parval_stem(14)
    
    emsoil   =parval_energy(11)

!   Initialize parameters and initial state ---------------------------------------
!   Initialize C parameters ----------------------------------------------
    LAI=1.0
 
    stor_use=Storage/times_storage_use
    accumulation=0.0
    Nscal_vcmax=1.0
    Vcmax_n=Vcmax*1.0e-6  ! change from umol/m2/s to mol/m2/s
    Jmax_n= 1.67*Vcmax_n

!   Initialize N parameters -----------------------------------------------
    NSN=1.0     ! Initial: 6.0
    NSNcrt1=0.10
    MineralN= 1.2
    NdemandT=0.0
    Ndeficit=0.0
    CN=CNini
    Npool=Cpool/CNini
    Ccost_Nacq=0.0
    Ndeficit=0.0
    
!   Initialize soil thermal dynamics -------------------------------------
    stemp=10.
    sftmp =10.
    Tsnow = -5.
    snow_depth=0.0
    water_free=0.0    
    dcount=50.
    zwt=0.0
    G=50.0
    ice_tw=0.0

!   Initialize soil_water dynamics -------------------------------------
    Rsoil=200
!   Initialize date ------------------------------------------------------
    m=1
    iyr=0
    idays=0
    day_of_year=365
    day_seq_spinup=0
    start_growing_season=91  ! 1st April
    end_growting_season=304   ! 31th October

    
    if (tec_diagn) print*,"tec_diag 1:"
    
!   Start of simulation for the whole time period  -------------------------
    do yr=1,year_length+yr_spinup  ! how many years
!        print*,"years:",year_length,yr_spinup,yr
!        stop
        
	
		
!   Assign number of days to a year ----------------------------------------
        iyr=iyr+1
        if(iyr>year_length)iyr=1
        if(MOD(first_year+iyr-1,4).eq.0)then
            day_of_year=366
            start_growing_season=92
            end_growting_season=305
        else
            day_of_year=365
            start_growing_season=91
            end_growting_season=304
        endif
!   Initialize yearly data --------------------------------------------------
        GDD5=0.0
        sene_count=0
        onset=0
        phenoset=0
        gpp_yr=0.0
        NPP_yr=0.0
        ANPP_yr=0.0
        
        NEE_yr=0.0
        Resp_eco_yr=0.0
        Resp_auto_yr=0.0
        Resp_hetero_yr=0.0
        NPP_L_yr=0.0
        NPP_R_yr=0.0
        LAI_yr=0.0
        
        transp_yr=0.0
        evap_yr=0.0
        runoff_yr=0.0
        water_leach_yr=0.0
        swater_slw_yr=0.0
        swater_slw_gs=0.0
        
        waterF_yr=0.0
        
!        print*,"stor_use 2:",stor_use,Storage,times_storage_use
        
        do days=1,day_of_year !the day of a year
            idays=idays+1
!            print*,"days:",days
            
!   Setting for N fertilization ----------------------------------------
!   Using ambient data to run equilibiurm, elevated only for the last cycle
            
!            if(yr>yr_spinup+1.and.(days==75.OR.days==105))then
!                MineralN=MineralN+Nfert     !(5.6 gN/yr/m2,N fertiliztion in March and Apr)
!            endif
            
!   Initialize some parameters ------------------------------------------
            tair_dmean=0.0   ! daily 
            rain_d=0.0     ! daily
            Resp_auto=0.0
            Resp_hetero=0.0
            Resp_eco=0.0
            NEE=0.0
            
            NEE_d=0.0
            GPP_d=0.0
            Resp_eco_d=0.0
            foliageC=0.0
            
            swc_d=0.0
            
!   Initialize C pools --------------------------------------------------

            NSNmin=0.1   ! NSCmin/CNini(1)*2.0

            
            if (do_snow) then 
                call snow_d(lat,yr,days,  &  ! input: space and time
                        &   rain_d,tair_dmean,  &   ! input: rain and temperature
                        &   fa,fsub,rho_snow,decay_m,   &   ! input: snow related parameters
                        &   dcount,snow_depth,melt)  ! outputs and/or inputs
            endif
      
            hour_of_day=24 !how many times a day,24 means every hour

            do i=1,hour_of_day
!                print*,"hour:",i,rcp
                
!   Set time line --------------------------------------------------------
!   Repeat forcing data for the whole time period
!   Once reach the end of forcing data, set data line (m) to be 1 to reuse the forcings. 
                if(m.gt.lines) m=1 
                if((yr.eq.(yr_spinup+1)).and.(days.eq.1).and.(i.eq.1)) then 
                    m=1
                else if ((yr.le.yr_spinup).and.(m.gt.lines)) then
                    m=1
                else
                endif
                
				!!! YZhou: add the following to add water to initial swc condition
				!if((yr.eq.(yr_spinup+1)).and.(days.eq.1).and.(i.eq.1)) then
				!	swc(7:10) = swc(7:10)*(1+0.6)  ! increase the water contents in deeper layers after spin up 
				!	print*, "Increase the water contents in deeper layers after spin up "
				!endif
				!!! YZhou: end of adding
				!!! YZhou: add following to recharge the deep soil 
				if (rechargeDS) then ! recharage the deep soil and keep it as a constant
					swc(10) = 0.3 
				endif
				!!! YZhou: end of adding
				
                year =year_seq(m)
!                print*,"size_year",size(year_seq)
!                stop
                doy  =doy_seq(m)
                hour =hour_seq(m)+1
                
!               print*,"Time:",year,doy,hour,hour_seq(m)

!   Prepare forcing data ------------------------------------------------- 
                if(rcp.eq.0)then
                    Tair=forcing_data(1,m)   ! Tair(degree)
                    co2ca=CO2treat*1.0E-6 
!                    print*,"Run under scenario 0."
                else if(rcp.eq.85)then
                    Tair=forcing_data(9,m)
                    co2ca=co2data(3,m)*1.0E-6
                    
!                    print*,"Run RCP8.5 scenario.",co2ca,Tair
                else if(rcp.eq.45)then
                    Tair=forcing_data(10,m)
                    co2ca=co2data(2,m)*1.0E-6
!                    print*,"Run RCP4.5 scenario."
                else if(rcp.eq.26)then
                    Tair=forcing_data(11,m)
                    co2ca=co2data(1,m)*1.0E-6
!                    print*,"Run RCP2.6 scenario."
                else
                    print*,"Warning: RCP scenario is incorrect!!!"
                    stop
                endif
                
                
                Tsoil=forcing_data(2,m)    ! Tsoil(degree)
                
                RH=forcing_data(3,m)        ! Relative humidity (%)
                Dair=forcing_data(4,m)       ! Saturation vapour pressure deficit (kPa)
                rain=forcing_data(5,m)    ! rainfall amount per hour (mm)
                wind=ABS(forcing_data(6,m))     ! wind speed m s-1
                radsol=forcing_data(7,m)        ! Solar radiation W/m2

                Tair=AMAX1(-55.0,AMIN1(55.0,Tair))   ! range [-50,50]
                Tsoil=AMAX1(-50.0,AMIN1(50.0,Tsoil)) ! range [-50,50]
                RH=AMAX1(0.01,AMIN1(99.99,RH))  ! range [0.01,99.99]
                Dair=AMAX1(0.01,AMIN1(20.0,Dair))  ! range [0.01,20]
                rain=AMAX1(0.00,AMIN1(200.0,rain)) ! range [0.0,200.0]
                wind=AMAX1(0.01,AMIN1(40.0,wind))    ! range [0.01,40.0]
                radsol=AMAX1(0.01,AMIN1(1500.0,radsol))    ! range [0.01,1500]
                
!               Recalculate Dair due to change in Tair
                Dair=(esat(Tair)-esat(Tair)*RH/100.)/1000
             
                
!   Add treatment effects ----------------------------------------------  
                if (yr .gt. 1)then
                    Tair = Tair+Ttreat
                    Tsoil = Tsoil+Ttreat
                    co2ca=CO2treat*1.0E-6 ! CO2 concentration,ppm-->1.0E-6
                endif
                TairK=Tair+273.2
                
!   Update parameters -------------------------------------------
                tair_dmean=tair_dmean+Tair/24.0  ! daily mean temperature for phenology and snow subroutines
                rain_d=rain_d+rain  ! daily rainfall for snow subroutine
                
                rhocp=cpair*Patm*AirMa/(Rconst*TairK)
!   rhocp: unit of Pa/oC^2, a product of air specific heat capacity (cpair) and air density
!                air density = air pressure (Pa)/(Temperature*Specific gas constant for dry air)
!                air density = Patm/(TairK
                
                
                Cmolar=Patm/(Rconst*TairK)
                H2OLv=H2oLv0-2.365e3*Tair  ! adjust H2OLv by air temp
                slope=(esat(Tair+0.01)-esat(Tair))/0.01   
                psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
                
                
!   Simulation starts following ---------------------------------------  
!                print*,"Vcmax_n:",Vcmax_n
                call canopy(m,yr,yr_spinup,doy,hour,LAI,sthick,layern,phenoset,frlen,& ! input: site settings
                &   wrt_tsoil,dosoilexp,&   ! switches of writing out
                &   lat,parval_stem,parval_energy,par_main,par_n,&  ! grouped variables
                &   swater_slw,swater_sld,swc,wsmax,wsmin,wsmin_fr,swc_fr,& ! input: water scaler
                &   radsol,co2ca,Tair,TairK,Dair,wind,rain,& ! input: climate
                &   Rconst,rhocp,Cmolar,H2OLv,psyc,slope, &   ! input: constants
                &   G,Esoil,Rsoil,&  ! input: energy fluxes
                &   Vcmax_n,Jmax_n,&    ! input: photosynthesis variables
                &   stemp,sftmp,raero,&    ! output: soil temperature
                &   snow_depth,Tsnow,zwt,ice,waterv,ice_tw,&  ! input: snow properties
                &   gpp,evap,transp,fbeam,&   ! outputs
                &   Eabsorp_leaf,RsoilabT,Rsoilab,Hsoil,Etransp,Hcan,gsc_l,Tleaf,Troot)    ! outputs
                
                if(do_soilphy)then
                    call soil_temp(Rsoilab,RsoilabT,Esoil,&  ! input
                    &   Tair,Cmolar,fbeam,LAI,sthick,layern,frlen,&  ! input
                    &   wsmax,wsmin,swater_sld,swc,&  ! input
                    &   parval_stem,Esoil_res,&  ! input
                    &   sigma,emsoil,Rsoil, &  ! input
                    &   Rconst,Rhocp,psyc,slope, & ! input: constants
                    &   cpair,Patm,AirMa,H2OMw,H2OLv0,H2OLv,raero,zwt,&  ! input
                    &   snow_depth,Tsnow, & ! input
                    &   ice,waterv,ice_tw,&  ! both input and output
                    &   wrt_tsoil,&   ! input
                    &   Hsoil,G,stemp,sftmp,Troot) ! outputs 
                else
                    
                    if(radsol.gt.10.0)then
                        G=-50
                    else
                        G=50
                    endif
                    sftmp=Tair
                    stemp=Tsoil
                    Troot=Tsoil
                    Hsoil=RsoilabT-Esoil-G
                endif
 
                if (tec_diagn) write(*,*) "tec_diag 3:",MineralN
              
                call soil_water(swc,wsmax,wsmin,watersat,waterres,condsat,wlama,potcof,n_pot,& ! input
                &   parval_swater,transp,evap,radsol,& 
                &   sthick,rdepth,frlen,layern,&
                &   water_free,rain,snow_depth,tair_dmean,&
                &   zwt,ice,waterv,melt, &    ! both input and output
                &   runoff,swater_slw,swater_sld,relsat,wsmin_fr,swc_fr,Rsoil,&   ! outputs
                &   waterF,wpotent,water_leach)   ! outputs  
                 
                if (tec_diagn) print*,"tec_diag 4:"

                call plant_respiration(vegetype,Tleaf,sftmp,Troot,Cpool,GPP,LAI,swater_slw,gcostpro,mresp20,   &   ! input
                &   Q10paccrate,    &
                &   RmainT,Rmain,RgrowthT,Rgrowth)  ! output

                call plant_growth(vegetype,GPP,RmainT,RgrowthT,Ccost_Nacq, &   ! input
                &   LAI,swater_slw, &   ! input
                &   NPPT,NPP,growthT,growth)    ! output
                
                call litter_fall(Cpool,TauC,swater_slw,Tair,Troot,   &   ! input
                &   litterfall,rootfall) ! output
                
                
                call hetero_resp(vegetype,outC,tauC,Growth,litterfall,rootfall,radsol,par_main,Q10,   & ! input
                &   swater_slw,stemp,basew4sresp,   &   ! input
                &   Cpool,  &   ! update
                &   Rh_pools)   ! output

                call Eco_N(Ndeposit,CNini,Cpool,NPP,outC,stemp,par_main,runoff,  &   ! input
                &   Npool,MineralN,NSN,NSNcrt1,Ndeficit,  &   ! update
                &   NdemandT,CN,NresorpT, &  ! output
                &   Nuptake,Nfix,Nvol,NimmobT,   &  ! output
                &   Nscal_vcmax,Nscal_autoresp,Ccost_Nacq)  ! output
                
                if (tec_diagn) print*,"tec_diag 7:",MineralN
                
!   Update variables ------------------------------------------------------
                LAI=Amax1(Cpool(1)/0.48*(SLA/10000),0.05) ! If LAI is 0.0, calculation based on LAI may have problems
!                print*,"LAI:",LAI,Cpool(1),SLA
                
                Resp_auto=RmainT+RgrowthT !+Ccost_Nacq
                stress_scal=Nscal_vcmax
                gdd_scal=(1+Amax1(Amin1((2500-GDD5)/2500,1.0),-1.0)*0.8)
                if(GDD5.lt.4000)then
                    LAI_scal=(2*(1.2-LAI))
                else
                    LAI_scal=1.0
                endif
                    
!                Vcmax_n= Vcmax*1.0e-6*stress_scal*gdd_scal*LAI_scal
                Vcmax_n= Vcmax*1.0e-6*stress_scal*gdd_scal
!                (1+Amax1(Amin1((2500-GDD5)/2500,1.0),-1.0)*0.8)  ! *stress_scal
!                print*,"Vcmax_n:",Vcmax_n,Vcmax,Nscal_vcmax
                
                Jmax_n= 1.67*Vcmax_n ! Weng 02/21/2011 Medlyn et al. 2002
!                water_free = water_free+rain !   Updated water available for soils
                
!   For screen showing, sum of the whole year data ----------------------------
                gpp_yr=gpp_yr+gpp
                NPP_yr=NPP_yr+NPPT
                ANPP_yr=ANPP_yr+NPP(1)+NPP(2)
                
                if(isnan(gpp))then
                    print*,'Warning: GPP is NAN!!!'
                    stop
                endif
                
!   For writing out ---------------------------------------------------------
!               Write out C pools --------------------
                if (yr.gt.yr_spinup) THEN
                    Resp_hetero=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)
                    Resp_eco=Resp_auto+Resp_hetero
                    NEE=Resp_eco-GPP
                    
                    NEE_d=NEE_d+NEE
                    GPP_d=GPP_d+GPP
                    Resp_eco_d=Resp_eco_d+Resp_eco
                    foliageC=foliageC+cpool(1)/24
                    
                    swc_d=swc_d+swc(1:5)/24
                    
                    if(wrt_cflux_yr.and.MCMC.eq.0)then
                        NEE_yr=NEE_yr+NEE
                        Resp_eco_yr=Resp_eco_yr+Resp_eco
                        Resp_auto_yr=Resp_auto_yr+Resp_auto
                        Resp_hetero_yr=Resp_hetero_yr+Resp_hetero
                        NPP_L_yr=NPP_L_yr+NPP(1)
                        NPP_R_yr=NPP_R_yr+NPP(3)
                        LAI_yr=LAI_yr+LAI/(day_of_year*24)
                    endif
              
!               Write out water fluxes -------------------
                    if(wrt_swater_yr.and.MCMC.eq.0)then
                        transp_yr=transp_yr+transp
                        evap_yr=evap_yr+evap
                        runoff_yr=runoff_yr+runoff
                        water_leach_yr=water_leach_yr+water_leach
                        swater_slw_yr=swater_slw_yr+swater_slw/(day_of_year*24)
                        if(days.ge.start_growing_season .and. days.le.end_growting_season)then  ! calculate the averaged soil water scaler of growing season (April to October)
                            swater_slw_gs=swater_slw_gs+swater_slw/((end_growting_season-start_growing_season)*24)  ! 274 is the length of growing season
                        endif
                        waterF_yr=waterF_yr+waterF
                    endif
                    
!               Write out soil temperature simulation ---
                    if (wrt_tsoil.and.MCMC.eq.0) then
                        write (3003,30031) m,sftmp,(stemp(j),j=1,10)
    !                    print*,"stemp:",stemp(1)
                    endif
30031               format(1(I8,","),10(f15.10,","),(f15.10))

!               Write out energy fluxes -------------------
!                    RsoilabT_yr=RsoilabT_yr+RsoilabT
!                    Eabsorp_leaf_yr=Eabsorp_leaf_yr+Eabsorp_leaf
!                    Esoil_yr=Esoil_yr+Esoil
!                    Hsoil_yr=Hsoil_yr+Hsoil
!                 
!                    Etransp_yr=Etransp_yr+Eabsorp_leaf


                    if (wrt_energy.and.MCMC.eq.0) then
                        write(3004,30041) m,RsoilabT,Eabsorp_leaf,Esoil,Hsoil,G,Etransp,Hcan,raero,gsc_l
                    endif
30041               format(1(I8,","),8(f15.10,","),(f15.10)) 

!               Write out N pools and fluxes ---------------
                    if (wrt_nflux.and.MCMC.eq.0) then
                        write(3005,30051) m,Nuptake,Nminer,Nimmob,NresorpT,Nfix,Nleach,Nrunoff,Nvol,NdemandT,Ndeposit, &
                                &   NewplantN(1),NewplantN(2),NewplantN(3),MineralN,NSN,&
                                &   Npool(1),Npool(2),Npool(3),Npool(4),Npool(5),Npool(6),Npool(7),Npool(8), &
                                &   Nscal_vcmax,Nscal_autoresp
!                    print*,"N leach:",Nleach,water_leach
                    endif
30051               format(1(I8,","),10(f15.12,","),14(f15.6,","),(f15.6))  

!               Write out soil water fluxes ---------
                    if (wrt_swaterflux.and.MCMC.eq.0) then
                        write (3007,30071) m,water_free,runoff,water_leach,waterF(1),waterF(2),waterF(3),waterF(4),&
                             &  waterF(5),waterF(6),waterF(7),waterF(8),waterF(9),wpotent(1),wpotent(2),wpotent(3),&
                             &  wpotent(4),wpotent(5),wpotent(6),wpotent(7),wpotent(8),wpotent(9),wpotent(10)
                    endif
30071               format(1(I8,","),21(f15.10,","),(f15.10))

!               Write out parameters for diagnose
                    if(wrt_pheno.and.MCMC.eq.0) then
                        write(3100,31001) m,phenoset,onset,sene_count,&
                            &   GDD5,add,store,storage,accumulation
                    endif
31001               format(4(I8,","),4(f15.10,","),(f15.10))

                endif
                
                
   
                m=m+1 ! move to next hour                
            enddo    ! end of hour_of_day -----------------------------------
            
!           write daily output 
            if(yr.gt.yr_spinup)then
                if(MCMC .ne.0 ) then
                    day_n=idays
                    simu_cflux_d(1,day_n)=NEE_d
                    simu_cflux_d(2,day_n)=GPP_d
                    simu_cflux_d(3,day_n)=Resp_eco_d 
                    simu_cflux_d(4,day_n)=FoliageC 
                    simu_swc_d(:,day_n)=swc_d
                endif
            
                if(MCMC.eq.0) then
                    if(wrt_cflux_d) then ! 
                        write(3001,30011) NEE_d,GPP_d,Resp_eco_d,FoliageC
                    endif
30011               format(3(f15.8,","),(f15.10))
                
                    if(wrt_swater_d) then
                        write(3002,30021) swc(1),swc(2),swc(3),swc(4),swc(5),swc(6),swc(7),swc(8),swc(9),swc(10) !! YZhou: addedd another five layers from 6 to 10
                    endif
30021               format (9(f15.10,","),(f15.10))

                endif

            endif

!   Update phenology for the next day -----------------------------------------
            if(tair_dmean.gt.5.0)GDD5=GDD5+tair_dmean
            
            if(tair_dmean.lt.5.0 .and. phenoset.eq.1 .and. days.gt.183) then
                sene_count=sene_count+1
            endif
!
!            
            if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
!                print*,"pheno"
                pheno=days
                phenoset=1
            endif
            
            if(sene_count.gt.5 .and. phenoset.eq.1) then
                phenoset=0
            endif
            
!            print*,"GDD5:",GDD5,onset,phenoset
            
        enddo      ! end of day_of_year

!   Updated storage and phenology for the next year ---------------------------
        storage=accumulation
        stor_use=Storage/times_storage_use
!        print*,"stor_use:",stor_use,Storage,times_storage_use
        accumulation=0.0
        onset=0
        
        if (yr.gt.yr_spinup) THEN  
!            Write out carbon fluxes --------------------------
            
            if(wrt_cflux_yr) then ! 
                write(3301,33011) iyr+first_year-1,GPP_yr,NPP_yr,NEE_yr,Resp_eco_yr,Resp_auto_yr, &
                    &   Resp_hetero_yr,NPP_L_yr,NPP_R_yr,LAI_yr
            endif
33011       format(1(I8,","),8(f15.6,","),(f15.6))

!           Write out water fluxes -------------------------
            if(wrt_swater_yr) then
                write(3302,33021) iyr+first_year-1,transp_yr,evap_yr,runoff_yr,water_leach_yr,swater_slw_yr,swater_slw_gs, &
                    &  waterF_yr(1),waterF_yr(2),waterF_yr(3),waterF_yr(4),waterF_yr(5) 
            endif
33021       format(1(I8,","),10(f15.10,","),(f15.10))


!           Write out energy fluxes --------------------------
!            if (wrt_energy) then
!                write(3304,33041) m,RsoilabT,Eabsorp_leaf,Esoil,Hsoil,G,Etransp,Hcan,raero,gsc_l
!            endif
!33041               format(1(I8,","),8(f15.10,","),(f15.10)) 




        endif
        
        
!   Screen showing ------------------------------------------------------------
        
        if (MCMC .NE. 1) then
            write(*,*)'annual_summary:',year,LAI,gpp_yr,NPP_yr,ANPP_yr
            if(NPP_yr.eq.0)then 
                print*,"Warning: yearly NPP is zero!!!"
!                stop
            endif           
        endif            

    enddo            !end of simulations multiple years
    
    if (tec_diagn) print*,"tec_diag 8:"
    
    return
end  
!   End of TECO_simu ===============================================================================
      
!****************************************************************************

!      a sub-model for calculating C flux and H2O flux of a canopy
!      adapted from a two-leaf canopy model developed by Wang Yingping
         
subroutine canopy(m,yr,yr_spinup,doy,hour,LAI,sthick,layern,phenoset,frlen,& ! input: site settings
            &   wrt_tsoil,dosoilexp,&   ! switches of writing out
            &   lat,parval_stem,parval_energy,par_main,par_n,&  ! grouped variables
            &   swater_slw,swater_sld,swc,wsmax,wsmin,wsmin_fr,swc_fr,& ! input: water scaler
            &   radsol,co2ca,Tair,TairK,Dair,wind,rain,& ! input: climate
            &   Rconst,rhocp,Cmolar,H2OLv,psyc,slope, &   ! input: constants
            &   G,Esoil,Rsoil,&  ! input: energy fluxes
            &   Vcmax_n,Jmax_n,&    ! input: photosynthesis variables
            &   stemp,sftmp,raero,&    ! output: soil temperature
            &   snow_depth,Tsnow,zwt,ice,waterv,ice_tw,&  ! input: snow properties
            &   gpp,evap,transp,fbeam,&   ! outputs
            &   Eabsorp_leaf,RsoilabT,Rsoilab,Hsoil,Etransp,Hcan,gsc_l,Tleaf,Troot)    ! outputs
       

    implicit none
    
    logical,parameter::can_diagn=.false.
    
!   inputs
    integer m,yr,yr_spinup,doy,hour,layern,phenoset,par_n
    real LAI,sthick(10),frlen(10) 
    logical wrt_tsoil,dosoilexp
    real par_main(par_n),parval_stem(20)
    real parval_energy(13)
    real swater_slw,swater_sld(10),swc(10),wsmax(10),wsmin(10),wsmin_fr,swc_fr
    real radsol,co2ca,Tair,TairK,Dair,wind,rain
    real G,Esoil,Esoil_res,Rsoil
    real Rconst,rhocp,Cmolar,psyc,slope
    real Vcmax_n,Jmax_n
    real stemp(10)
    real snow_depth,Tsnow,zwt,ice(10),waterv(10),ice_tw
    real Vcmax_par(3)

!   outputs
    real gpp,evap,transp,Eabsorp_leaf,RsoilabT,Hcan,gsc_l,Tleaf,Troot
    
!   internal variables
    real Acanop,alpha,airMa,coszen,coszen_cal,cpair,conKo0,chi
    real Dheat,Ds0,Etransp,Ekc,Eko,esat,eairP,emleaf,emsoil
    real evap_pot,extkbm(3),extkdm(3),extkU,extKb,extKd,fbeam
    real Gaussx(5),Gaussw(5),Gaussw_cum(5),gam0,gam1,gam2,gsc0,gddonset
    real Hsoil,H2OLv0,H2OLv,H2OMw
    real lat,Kc_ref
    real Patm,pi,Qcan(3,2),o2ci,QLair,QLleaf
    real Rsoilab(3),raero
    real rhoL(3),rhoS(3),reffbm(3),reffdf(3),Radabv(2)
    real sigma,sftmp,swater_crt1,swater_crt2,swater_crt3,swater_crt4,stom_n
    real tauL(3),tew,rew,transp_pot,Temp_ref
    real wleaf,xfang
    
    
    integer jrain,i,j,k

    character*80 commts
!     Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!     5-point
    data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
    data Gaussw/0.118430025,0.2393144,0.2844444,0.2393144,0.1184635/
    data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/

    
!   Constants --------------------------------
    sigma=5.67E-08  ! Steffan Boltzman constant; unit: W/m2/K4
    Temp_ref=293.2  ! Temperature reference; 25OC
!    H2OLv0=2.50E+06 ! Latent heat of H2O; unit: J/kg
    airMa=2.90E-02  ! Molar mass of air; unit: kg/mol
    H2OMw=1.80E-02  ! Molar mass of H2O; unit: kg/mol
    Dheat=2.15E-05  ! Molecular diffusivity for heat; unit: m2/s
    Patm=101325     ! Atmospheric pressure; unit: Pa
    cpair=1010      ! Heat capacity of air; unit of J/kg/K

    stom_n=par_main(3)
    Ds0=par_main(4)
    gddonset=par_main(14)
    Vcmax_par(1)=par_main(23)
    Vcmax_par(2)=par_main(24)
    Vcmax_par(3)=par_main(25)

    tauL(1)  =parval_energy(1)
    rhoL(1)  =parval_energy(2)
    rhoS(1)  =parval_energy(3)
    tauL(2)  =parval_energy(4)
    rhoL(2)  =parval_energy(5)
    rhoS(2)  =parval_energy(6)
    tauL(3)  =parval_energy(7)
    rhoL(3)  =parval_energy(8)
    rhoS(3)  =parval_energy(9)
    emleaf   =parval_energy(10)
    emsoil   =parval_energy(11)
    wleaf    =parval_energy(12)
    gsc0     =parval_energy(13)
    
    Esoil_res=0.0
        
!   Calculated for grassland, according to Box 4 at the FAO website: http://www.fao.org/docrep/x0490e/x0490e06.htm#chapter%202%20%20%20fao%20penman%20monteith%20equation         
    raero=Amin1(200/wind,2000.0)   ! aerodynamic resistance                        


    if(can_diagn) print*,"can_diag 1:"
    
!   Calculate beam fraction in incoming solar radiation
    call  fbeam_cal(doy,hour,lat,radsol,fbeam)
!   Calculate cos zenith angle of sun
    coszen=coszen_cal(doy,lat,hour) 

!   Calculate soil albedo for NIR as a function of soil water (Garratt pp292)
    if(swater_sld(1).gt.0.5) then
        rhoS(2)=0.18
    else
        rhoS(2)=0.52-0.68*swater_sld(1)   ! use water scaler
    endif

    
    eairP=esat(Tair)-Dair*1000                !air water vapour pressure, unit of Pa
    if (eairP.le.0) then
        eairP=1.0
!        print*,"Warning: unreasonable water vapour deficit!!!"
    endif

    radabv(1)=0.5*radsol                 !(1) - VIS radiation/ PAR
    radabv(2)=0.5*radsol                 !(2) - NIR
 
    if(can_diagn) print*,"can_diag 2:"
!   print*,"gpp0:",gpp,Acanop
!   Call multilayer model of Leuning - uses Gaussian integration but radiation scheme   
!   print*,"Vcmax_n:",Vcmax_n
    call xlayers(m,yr,yr_spinup,phenoset,  & ! input
            &  radabv,fbeam,co2ca,Dair,eairP,Tair,TairK,wind,sftmp,Vcmax_par,& ! input     ! 
            &  wsmax,wsmin,swater_slw,swater_sld,LAI, & ! input
            &  coszen,alpha,Rconst,rhocp,Cmolar,psyc,slope,& ! input
            &  tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,& ! input
            &  sigma,emleaf,emsoil,Ds0,& ! input
            &  cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat,& ! input
            &  gsc0,stom_n,raero,Rsoil,& ! input
            &  Vcmax_n,Jmax_n,Kc_ref,conKo0,Ekc,Eko,o2ci,& ! input
            &  gam0,gam1,gam2,& ! input
            &  G, & ! input
            &  RsoilabT,Rsoilab,&  ! output
            &  Acanop,Etransp,QLleaf,QLair,extKb,&  ! output
            &  Esoil,Eabsorp_leaf,Hcan,gsc_l,Tleaf)  ! output
        
    if(can_diagn) print*,"can_diag 3:"
    
!   Calculate hourly gpp; 
!   Unit of Acanop: mol C/m2/s
!   Unit of gpp: g C/m2/h
    gpp=Acanop*3600.0*12.0

!   Potential transpiration; Unit of Etransp is W/m2 or J/m2/s; Unit of transp is mm/h
    transp_pot=AMAX1(Etransp/(2.501-0.00236*Tair)*3600.0/1.0e6,0.) ! mm H2O /hour
!   Actual transpiration, according to Novák, V., & Havrila, J. (2006).
    swater_crt2=0.67*wsmin_fr
    swater_crt1=swater_crt2+1/(-2.27*(transp*24)+17.5)  ! transp*24 for daily transp (mm/day)
    
    if (swc_fr.lt.swater_crt2) then
        transp=0.0
    else if (swc_fr.lt.swater_crt1) then
        transp=transp_pot*(swc_fr-swater_crt2)/(swater_crt1-swater_crt2)
    else
        transp=transp_pot
    endif
    
    if(transp_pot.ne.0) Esoil_res=Esoil_res+Etransp*(transp_pot-transp)/transp_pot
 

!    print*,swater_crt1,swater_crt2,swc_fr,transp,transp_pot,transp/transp_pot*100
    
!   Unit of Esoil is W/m2 or J/m2/s; Unit of evap is mm/h; similar units for transp
!   1 mm/day = 2.45 MJ/m2/day = 0.102 MJ/m2/h ~ 0.0417 mm/h
!   (2.501-0.00236*Tair) converted MJ/m/h to mm/h, meanwhile calibrated with air temperature (oC) 
!   Potential evaporation
    evap_pot=AMAX1(Esoil/(2.501-0.00236*Tair) *3600.0/1.0e6,0.)        
!   Actual evaporation, according to FAO: http://www.fao.org/docrep/x0490e/x0490e0c.htm#TopOfPage
!   rew: readily evaporable water, which is the maximum depth of water 
!        that can be evaporated from the topsoil layer without restriction (mm)
!   rew ranges between 0.04-0.10 for loamy sand and sandy loam soils, and between 0.02 and 0.12 for all soil types
!   tew: maximum accumulative depth of evaporation (depletion) from the soil surface layer (mm)
!   tew ranges between 0.09-0.20 for loamy sand and sandy loam soils and between 0.06-0.29 for all soils 
    if(dosoilexp)then
        rew=(wsmax(1)-wsmin(1))*0.6
    else
        rew=0.02
    endif
    tew=wsmax(1)-0.5*wsmin(1) ! initial 0.5   ! 0.30
    swater_crt3=wsmax(1)-rew    ! 0.35
    swater_crt4=wsmax(1)-tew        ! 0.15
    if (swc(1).gt.swater_crt3) then
        evap=evap_pot
    else if (swc(1).gt.swater_crt4) then
        evap=evap_pot*(swc(1)-swater_crt4)/(swater_crt3-swater_crt4)
    else
        evap=0
    endif
        
    if(evap_pot.ne.0)  Esoil_res=Esoil_res+Esoil*(evap_pot-evap)/evap_pot


!    print*,"evap:",evap


    return
end  ! End of subroutine canopy
!   *******************************************************************************

!*******************************************************************
!     subroutine for soil moisture
subroutine soil_water(swc,wsmax,wsmin,watersat,waterres,condsat,wlama,potcof,n_pot,& ! input
                &   parval_swater,transp,evap,radsol,& 
                &   sthick,rdepth,frlen,layern,&
                &   water_free,rain,snow_depth,tair_dmean,&
                &   zwt,ice,waterv,melt, &    ! both input and output
                &   runoff,swater_slw,swater_sld,relsat,wsmin_fr,swc_fr,Rsoil,&   ! outputs
                &   waterF,wpotent,water_leach)   ! outputs 
!   All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
!   Notes: unit of waterv,swc,and ice are all v/v.
                
    implicit none
    
!   Switches
    logical,parameter:: do_hydr_redist=.false.
    
!   Inputs
    integer layern
    real swc(10),wsmax(10),wsmin(10),watersat(10),waterres(10),condsat(10)
    real parval_swater(6),evap,transp,radsol
    real sthick(10),rdepth,depth(10),frlen(10),bmroot
    real water_free,snow_depth,zwt,ice(10),waterv(10)
    real melt,rain,tair_dmean,alpha
   
!   Outputs
    real runoff,swater_slw,swater_sld(10),wsmin_fr,swc_fr,Rsoil
    real waterF(10),WPotent(10),water_leach
 
!   Internal variables 
    integer nfr,i,j,k,i_transp,count_transp
    real Crt,Crt_ad,c(10),Fc(10),Wtheta(10),condval(10),condb(10),Dtran  !! YZhou: added the Fc and Wtheta (Ryel et al 2002)
    real Frlenx(10),HRP_tem
    real mpara(10),n_pot(10),plantup(10),potcof(10)
    real relsat(10),relwater(10)
    real Tr_ratio(10),tr_allo
    real wsc(10),wsc_emp(10),wsc_emp_acc(10),wlama(10) ! Hydraulic redistribution
    real WPotentDT, infirate
    real vtot
    real az,phi,thetasmin,zmax,zthetasmin,zwt1,zwt2,zwt3
 
    
!    Crt=parval_swater(1)
!    wlama=parval_swater(2)
!    potcof=parval_swater(3)
!    condb=parval_swater(4)
!    n_pot=parval_swater(5)

    count_transp=0
    DEPTH=0.0
    relsat=-9999.0
    relwater=-9999.0
    waterF=-9999.0
    swater_sld=-9999.0
    WPotent=-9999.0
    tr_ratio=-9999.0
    plantup=-9999.0

!   Determine which layers are reached by the root system. 
!   Layer volume (cm3)
    DEPTH(1)=sthick(1)
    do i=2,layern
        DEPTH(i)=DEPTH(i-1)+sthick(i)
    enddo
    
    do i=1,layern
        IF(rdepth.GT.DEPTH(i)) nfr=i+1
    enddo

    IF (nfr.GT.10) nfr=10

!   Melting and rainfall---------------------------------------------------------------------

    if (tair_dmean .lt. -4.) then
        water_free=water_free+melt/24
    else
        water_free=water_free+melt/24+rain
    endif

    if (ice(1) .gt. 0.0) then
        !water_free = 0.0
    endif
	!print*, water_free
    print*, "reach here..., and swc(10) is", swc(10)
!   Rainfall (mm h–1) was added to dWi/dt for the uppermost layer only, Ryel et al (2002)
    ! calculation the infliation rate based on Eq. 8.10 in Elements of Physical Hydrology
	if (water_free .gt. 0) then
		n_pot(1) = 2.06
		mpara(1)=1-1/n_pot(1) 
		Wtheta(1) = swc(1)
		Wpotent(1)= -((((watersat(1)-waterres(1))/(Wtheta(1)-waterres(1)))**(1/mpara(1)))-1)**(1/n_pot(1))/alpha 
		print*, "WP = ", Wpotent(1), " ; ws -wr =", watersat(1)-waterres(1), "swc = ", Wtheta(1), "swc- wr =",Wtheta(1)-waterres(1)
		print*,"m=",mpara(1), (((watersat(1)-waterres(1))/(Wtheta(1)-waterres(1)))**(1/mpara(1)))
		infirate = -condsat(1)*((-Wpotent(1)-0)/(sthick(1))+1) ! infiltration rate cm h-1
		print*, "infi rate = ", infirate, " ;top soil SWC and WP are", swc(1), "and", Wpotent(1)
		if (infirate*1 .lt. water_free) then
			swc(1)= swc(1)+infirate/10/sthick(1) !cm3/cm3 !water_free comes from rainfall and snowmelt 
			runoff=water_free-infirate*1 ! 1 represents one hour step
			print*, "Runoff = ", runoff, " ;infi rate = ", infirate, " ;top soil SWC and WP are", swc(1), "and", Wpotent(1)
		else
			print*,"adding water ",water_free/10/sthick(1), " to top soil from ", swc(1)
			swc(1)= swc(1)+water_free/10/sthick(1)			
		endif
	endif
    !   Unsaturated soil water flow ------------------------------------- 
!   Calculation of soil water potential from soil water content, according to Ryel et al. (2002)
!   WPotent=((swc-waterres)/(watersat-waterres)**((1/-wlama)/potcof)     
!   condval (unsaturated soil hydraulic conductivity), is calculated using van Genuchten (1980)
!   wpotent calculated according to Brooks and Corey (1964) and Moreno‐de las Heras et al. (2016)
!   wpotent may also be calculated according to van Genuchten (1980) and Rse (1997)
!   n_pot = 2.0 ! decrease n_pot there will be a smaller decrease in swc when swc is high
    do i=2,layern
		!relsat(i)=Amax1(0.0001,Amin1(0.9999,(swc(i)-waterres(i))/(watersat(i)-waterres(i)))) ! relative soil water saturation, dimensionless, Eq(4) in Ryel et al (2002)
        !relsat(i+1)=Amax1(0.0001,Amin1(0.9999,(swc(i+1)-waterres(i+1))/(watersat(i+1)-waterres(i+1)))) ! relative soil water saturation, dimensionless
        n_pot(i) = 2.0 !! YZhou: paramters for type silt loam GE3 from van Genuchten (1980) !! need to be put at the begining or somewhere else
		n_pot(i-1) = 2.0
		alpha = 0.005
		mpara(i)=1-1/n_pot(i) !!YZhou: change to the original value in Ryel et al (2002), fitting parameter for soil water retention curve		!n_pot(i) 
        mpara(i-1)=1-1/n_pot(i-1) !n_pot(i-1) 
        ! soil water matric potential (MPa; 1Mpa = 10200 cm)
		Wtheta = swc ! could be removed and use swc directly
		Wpotent(i) = -1/alpha * (((watersat(i) - waterres(i))/(Wtheta(i) - waterres(i)))**(1/mpara(i))-1)**(1/n_pot(i))
		Wpotent(i-1) = -1/alpha * (((watersat(i-1) - waterres(i-1))/(Wtheta(i-1) - waterres(i-1)))**(1/mpara(i-1))-1)**(1/n_pot(i-1))

		!Wpotent(i)= -((((watersat(i)-waterres(i))/(Wtheta(i)-waterres(i)))**(1/mpara(i)))-1)**(1/n_pot(i))/alpha 
		!Wpotent(i-1)= -((((watersat(i-1)-waterres(i-1))/(Wtheta(i-1)-waterres(i-1)))**(1/mpara(i-1)))-1)**(1/n_pot(i-1))/alpha 

		! relative saturation
		relsat(i) = (Wtheta(i)-waterres(i))/(watersat(i)-waterres(i))
		relsat(i-1) = (Wtheta(i-1)-waterres(i-1))/(watersat(i-1)-waterres(i-1))
		! YZhou: Based on Eq(3) in Ryel et al (2002), the unsaturated soil hydraulic conductivity was expressed as a function of saturated hydraulic conductivity Ks and volumetric water content using van Genuchten (1980).
		! condsat: Ks in Eq(3), saturated hydraulic conductivity
		! relsat: Si in Eq(3), relative saturation
		! unsaturated soil hydraulic conductivity
        condval(i)=condsat(i)*relsat(i)**0.5 *(1-(1-relsat(i)**(1/mpara(i)))**mpara(i))**2  ! cm h-1     
		!print*,"cond = ",condval(i), relsat(i), Wtheta(i)
        !!YZhou: condval(i-1)=condsat(i-1)*relsat(i-1)**0.5 *(1-(1-relsat(i-1)**(1/mpara(i-1)))**mpara(i-1))**2 
        !!!WpotentDT =(WPotent(i)-WPotent(i-1)) ! cm
        !if (WpotentDT.ge.0) then
            !waterF(i)=Amax1(Amin1(condval(i),condval(i)*WpotentDT,(swc(i)-waterres(i))*(sthick(i)*10), &
                    !&   (watersat(i-1)-swc(i-1))*(sthick(i-1)*10)),0.0) ! unit: cm/h. Since WPotent is negative, waterF is negative.
		!!!!waterF(i) = condval(i)*(WpotentDT/((sthick(i)+sthick(i-1))/2)+1)  !cm h -1 = cm h-1 * (cm / cm)
        !else
            !waterF(i)=Amin1(Amax1(-condval(i-1),condval(i-1)*WpotentDT,(swc(i)-watersat(i))*(sthick(i)*10),  &
                    !&   (waterres(i-1)-swc(i-1))*(sthick(i-1)*10)),0.0) ! unit: cm/h. Since WPotent is negative, waterF is negative.            
        !endif
		WpotentDT=(WPotent(i)-WPotent(i+1))/((sthick(i)+sthick(i+1))/2)+1
        if (WpotentDT.ge.0) then
            waterF(i)=Amax1(Amin1(condval(i),condval(i)*WpotentDT,(swc(i)-waterres(i))*(sthick(i)*10), &
                    &   (watersat(i+1)-swc(i+1))*(sthick(i+1)*10)),0.0) ! unit: cm/h. Since WPotent is negative, waterF is negative.
        else
            waterF(i)=Amin1(Amax1(-condval(i+1),condval(i+1)*WpotentDT,(swc(i)-watersat(i))*(sthick(i)*10),  &
                    &   (waterres(i+1)-swc(i+1))*(sthick(i+1)*10)),0.0) ! unit: cm/h. Since WPotent is negative, waterF is negative.            
        endif
             
        swc(i)=swc(i)+waterF(i)/(sthick(i)) !cm3 cm-3, update swc for each hourly step
        swc(i-1)=swc(i-1)-waterF(i)/(sthick(i-1))  !cm3 cm-3
        print*,i,swc(i-1), swc(i), Wpotent(i-1),Wpotent(i),WpotentDT,condval(i),waterF(i),sthick(i)
		if (i .eq. 9) then
			print*,i,Wpotent(i), alpha, watersat(i),waterres(i),Wtheta(i),mpara(i),n_pot(i)
			print*,i-1,Wpotent(i-1), alpha, watersat(i-1),waterres(i-1),Wtheta(i-1),mpara(i-1),n_pot(i-1)
			print*,((watersat(i-1) - waterres(i-1))/Wtheta(i-1)- waterres(i-1))
			
		endif
        if(swc(i).lt.waterres(i)) then
            print*,i,swc(i),waterF(i),WPotent(i),WPotent(i-1),WPotentDT
            print*,"Warning: Something is wrong with unsaturated soil water conductance!!!"
            stop
        endif
    enddo
    
    

!   Runoff -------------- !!!!YZhou: may need to put to the end of soil moisture module after infliation 
    !runoff=water_free * 0.001
    !water_free=water_free-runoff
	!if (runoff .gt. 0) then
	!	print*,"runoff:",runoff,wsc_emp(1), water_free,wsmax(1)-swc(1),sthick(1)
	!endif
    
!   Evaporation, depletion of surface soil water content ----------------------------
!   evap, unit of mm/h
!    print*,"swc(1):",swc(1),swc(1) - evap/(sthick(1)*10.0),waterres(1)
    
    do i=1,layern
        relwater(i)=Amax1(0.00,(swc(i)-waterres(i))*sthick(i)*10.0) ! unit: mm
    enddo
    
!   Assume evap can depelte water at most in the first two layers of soils
 
    if (evap.le.relwater(1)) then
        swc(1)=swc(1)-evap/(sthick(1)*10.0)
    else if(evap.lt.(relwater(1)+relwater(2))) then
!        print*,"evap:",evap,relwater(1)
        swc(2)=swc(2)-(evap-relwater(1))/(sthick(2)*10.0)
        swc(1)=waterres(1)
    else
        swc(2)=waterres(2)
        swc(1)=waterres(1)
!        print*,"swc(1):",swc(1),swc(2),relwater(1),relwater(2),evap
!        print*,"Warning: Surface soil water is not enough for evaporation!!!"
!        stop
    endif
    
!   Transpiration ------------------------------------------------------------
!   Redistribute transpiration to soil layers according to root biomass
!   Unit of transp: mm/h
    do i=1,layern
        relwater(i)=Amax1(0.00,(swc(i)-waterres(i))*sthick(i)*10.0) ! unit: mm
    enddo 
    
    tr_allo=0.0
    do i=1,nfr
        tr_ratio(i)=FRLEN(i)*relwater(i) 
        tr_allo=tr_allo+tr_ratio(i)
    enddo
    
!    if( (tr_allo.eq.0.0).or.(transp.gt.tr_allo)) then
!        print*,"relative water:",tr_allo,transp,relwater
!        print*,"Warning: not enough soil water for transpiration!!!"
!        stop
!    endif
   
    do i=1,nfr
        plantup(i)=transp*tr_ratio(i)/tr_allo
!        if(((swc(i)-waterres(i))*sthick(i)*10.0).lt.plantup(i)) then
!            print*,"swc:",i,swc(i),tr_ratio(i),tr_allo,plantup(i),waterres(i),relwater(i) 
!            print*,"Warning: soil water is not enough for transpiration!!!"
!            stop
!        endif
        swc(i)=AMAX1(swc(i)-plantup(i)/(sthick(i)*10.0),waterres(i) )
    enddo
    
    
!   Hydraulic redistribution  ------------------------------------------------------------------
!   According to Ryel et al. (2002), Oecologia, 130(2), 173-184.
!   Equation: Hi = Crt*Sum(WPotent(j)-WPotent(i))*max(c(i),c(j))*(R(i)*R(j)/(1-Rx))*Dtran
!   Hi: water move into soil layer i from all other layers 

    if(do_hydr_redist) then
        if(radsol.gt.10)then ! radsol > 10 indicates day time
            Dtran=0  ! YZhou: changed 0 and 1 to indicate HR water movement in nighttime
        else
            Dtran=1  
        endif
        Crt_ad = 0.097! YZhou: orginal was Crt
        do i=1,10
            do j=1,10
				Wtheta(i) = (watersat(i)-waterres(i))/((1+abs(0.2828 * Wpotent(i))**1.40)**(1-1/1.40))
				Wtheta(j) = (watersat(j)-waterres(j))/((1+abs(0.2828 * Wpotent(j))**1.40)**(1-1/1.40))
                if(Wtheta(i).gt.Wtheta(j))then
                    Frlenx(i)=Frlen(i)
                else
                    Frlenx(i)=Frlen(j)
                endif
				! YZHou: adding c calculation, replacing c with Fc in case there is conflict
				Fc(i) = 1/(1+(WPotent(i)/(-1))**3.22)
				Fc(j) = 1/(1+(WPotent(j)/(-1))**3.22)
				! YZhou: end of adding
                HRP_tem = Crt_ad*(WPotent(j)-WPotent(i))*Amax1(Fc(i),Fc(j))*(Frlen(i)*Frlen(j))/(1-Frlenx(i))*Dtran
				!print*, "HRP_tem:", Crt_ad, HRP_tem, "WP and c of ", j, " are ",WPotent(j),  Fc(j), ";of ",i, " are ",WPotent(i),  Fc(i)
                swc(j)=swc(j)-HRP_tem/sthick(j)
                swc(i)=swc(i)+HRP_tem/sthick(i)
!				print*, "HR:", Dtran, HRP_tem, "from ", j, " Fr is ",Frlen(j),  swc(j), "to ", i," Fr is ",Frlen(i), swc(i)
				
            enddo
        enddo
    endif
    

!   Leaching ----------------------------------------------------------------
    water_leach=Amax1((swc(layern)-wsmax(layern))*sthick(layern)*10*0.0005,0.0)  ! unit of mm/h 
	!!YZhou: change the wsmax (water holding capacity), definitely not wsmin!!

    swc(layern)=swc(layern)-water_leach/(sthick(layern)*10)
!    print*,swc(layern),water_leach
    
!    print*,layern,swc(layern),(sthick(layern))
!   Soil water scalers for other subroutines --------------------------------
    swater_slw=0.0
    wsmin_fr=0.00
    swc_fr=0.0
    do i=1,layern       
        swater_sld(i)=(swc(i)-wsmin(i))/(wsmax(i)-wsmin(i))
        swater_sld(i)=AMIN1(1.0,AMAX1(0.0,swater_sld(i)))   ! linear (l) water scaler (s) of each layer of soil depth (d)
        swater_slw=swater_slw+swater_sld(i)*frlen(i)  ! linear (l) soil water scaler (s) weighted by the whole (w) root system 
        wsmin_fr=wsmin_fr+wsmin(i)*frlen(i)   ! root biomass weighted wilting point
        swc_fr=swc_fr+swc(i)*frlen(i)  ! root biomass weighted soil water content
    enddo
    
!    rsoil=30.*exp(0.2/AMIN1(AMAX1(1-swater_sld(1),0.20),1.0))   ! initial: 0.2, EH mod
    rsoil=500
    
!   Water table --------------------------------------------------------------

    vtot = ((waterv(1)+ice(1))*sthick(1)+(waterv(2)+ice(2))*sthick(2)+ &
        &   (waterv(3)+ice(3)*sthick(3)))*0.01+water_free

    
    phi = 0.56   !soil porosity   mm3/mm3   the same unit with theta
    zmax = 300   !maximum water table depth   mm
    thetasmin = 0.25    !minimum volumetric water content at the soil surface   cm3/cm3
    zthetasmin = 100     !maximum depth where evaporation influences soil moisture   mm
    az = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1

    zwt1 = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
    zwt2 = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
    zwt3 = vtot-phi*zmax                                   
    if ((zwt1 .ge. -100) .and. (zwt1 .le. 0))   zwt=zwt1  !the non-linear part of the water table changing line
    if (zwt2 .lt. -100)                         zwt=zwt2  !the linear part of the water table changing line
    if (phi*zmax .lt. vtot)                     zwt=zwt3  !the linear part when the water table is above the soil surface 


    return
end  !   End of soil_water subroutine
!   End of soil_water subroutine ==============================================================================  

!   plant growth subtroutine ==================================================================================
subroutine plant_growth(vegetype,GPP,RmainT,RgrowthT,Ccost_Nacq, &   ! input
            &   LAI,swater_slw, &   ! input
            &   NPPT,NPP,growthT,growth)    ! output
  
!   Define variables --------------------------------------------------------
    implicit none
!   input
    integer vegetype
    real GPP,RmainT,RgrowthT,Ccost_Nacq
    real LAI,swater_slw
!   output
    real NPPT,NPP(3),growthT,growth(3)
    
!   internal
    real growth_prop(3),growth_prop_base(3),weight_allo,growth_scalL

!    print*,"Ccost_Nacq:",Ccost_Nacq
 
!   Plant biomass growth ---------------------------------------------------  
!   According to CTEM model, Arora, V. K., & Boer, G. J. (2005) and Arora, V. K. (2003).
    GrowthT=GPP-RmainT-RgrowthT-Ccost_Nacq ! GrowthT can be negative, means C cost is higher than C gain.
    growth_scalL=Amax1(1-LAI/2.0,0.0)  ! Scaler of light availability; initial: Amax1(1-LAI/4.5,0.0)
    weight_allo=1.0  ! 
    
    if(vegetype.eq.1)then   ! 1 for grassland, without stem
        growth_prop_base(1)=0.5   ! initial: 0.45
        growth_prop_base(2)=0.0
        growth_prop_base(3)=1-growth_prop_base(1)-growth_prop_base(2)
        
        growth_prop(1)=(growth_prop_base(1)+weight_allo*growth_scalL)/(1+weight_allo*(1+growth_scalL-swater_slw))
        growth_prop(2)=growth_prop_base(2)
        growth_prop(3)=(growth_prop_base(3)+weight_allo*(1-swater_slw))/(1+weight_allo*(1+growth_scalL-swater_slw))

    else       
        growth_prop_base(1)=0.40
        growth_prop_base(2)=0.05
        growth_prop_base(3)=1-growth_prop_base(1)-growth_prop_base(2)
        
        growth_prop(1)=growth_prop_base(1)/(1+weight_allo*(2-growth_scalL-swater_slw))
        growth_prop(2)=(growth_prop_base(2)+weight_allo*(1-swater_slw))/(1+weight_allo*(2-growth_scalL-swater_slw))
        growth_prop(3)=1-growth_prop(1)-growth_prop(2)
    endif
    
    Growth(1)=GrowthT*growth_prop(1)
    Growth(2)=GrowthT*growth_prop(2)
    Growth(3)=GrowthT*growth_prop(3)
    
!    print*,"Growth(1):",Growth(1),GrowthT,growth_prop(1),GPP,RmainT,RgrowthT,Ccost_Nacq

!   Total NPP and NPP allocation --------------------------------------------------------------
!   NPP allocation
    NPP(1)=Growth(1) !  Amax1(Growth(1),0.0)
    NPP(2)=Growth(2) !  Amax1(Growth(2),0.0)
    NPP(3)=Growth(3) !  Amax1(Growth(3),0.0)
    
    NPPT=sum(NPP)
    
!    print*,"NPP:",NPPT,GrowthT,GPP,RmainT,RgrowthT,Ccost_Nacq
    
    return
end  
!   End of plant growth subroutine  ========================================================

!   ========================================================================================
subroutine litter_fall(Cpool,TauC,swater_slw,Tair,Troot,   &   ! input
                &   litterfall,rootfall) ! output
    implicit none
!   input
    real Cpool(8),TauC(8),swater_slw,Tair,Troot
!   output
    real litterfall,rootfall
!   internal
    real litterfall_scal,LFall_scalW,LFall_scalWmax,LFall_scalWb
    real LFall_scalT,LFall_scalTmax,LFall_scalTcrt,LFall_scalTb,Tcold
    real rootfall_scal,RFall_scalW,RFall_scalWmax,RFall_scalWb
    real RFall_scalT,RFall_scalTmax,RFall_scalTcrt,RFall_scalTb
    
!   Litterfall, according to CTEM model, Arora, V. K., & Boer, G. J. (2005)
    LFall_scalTmax=0.15/24
    LFall_scalTb=3.0
    Tcold=5.0  ! initial 5.0
    if(Tair.ge.Tcold)then
        LFall_scalTcrt=1.0
    else if (Tair.gt.(Tcold-5.0)) then
        LFall_scalTcrt=(Tair-(Tcold-5.0))/5.0
    else
        LFall_scalTcrt=0.0
    endif
    LFall_scalT=LFall_scalTmax*(1-LFall_scalTcrt)**LFall_scalTb
    
    LFall_scalWmax=0.005/24   ! initial: 0.005/24
    LFall_scalWb=3.0  ! initial: 3.0
    LFall_scalW=LFall_scalWmax*(1-swater_slw)**LFall_scalWb
    
    litterfall_scal=1.0/tauC(1)+LFall_scalW+LFall_scalT    
      
    litterfall=Cpool(1)*litterfall_scal

!   Root turnover according to litterfall
    RFall_scalTmax=0.15/24
    RFall_scalTb=3.0
    Tcold=5.0  ! initial 5.0
    if(Troot.ge.Tcold)then
        RFall_scalTcrt=1.0
    else if (Troot.gt.(Tcold-5.0)) then
        RFall_scalTcrt=(Troot-(Tcold-5.0))/5.0
    else
        RFall_scalTcrt=0.0
    endif
    RFall_scalT=RFall_scalTmax*(1-RFall_scalTcrt)**RFall_scalTb
    
    RFall_scalWmax=0.005/24   ! initial: 0.005/24
    RFall_scalWb=3.0  ! initial: 3.0
    RFall_scalW=RFall_scalWmax*(1-swater_slw)**RFall_scalWb
    
    rootfall_scal=1.0/tauC(3)+RFall_scalW+RFall_scalT    
      
    rootfall=Cpool(3)*rootfall_scal

    
    return
end ! end of litter_fall subroutine
!   ========================================================================================

!   ========================================================================================
!     carbon transfer according to Xu et al. 2007 
subroutine Eco_N(Ndeposit,CNini,Cpool,NPP,outC,stemp,par_main,runoff,  &   ! input
            &   Npool,MineralN,NSN,NSNcrt1,Ndeficit,  &   ! update
            &   NdemandT,CN,NresorpT, &  ! output
            &   Nuptake,Nfix,Nvol,NimmobT,   &  ! output
            &   Nscal_vcmax,Nscal_autoresp,Ccost_Nacq)  ! output
    
!   Define variables
    implicit none
!   input
    real Ndeposit,CNini(8),Cpool(8),NPP(3),outC(8)
    real stemp(10),par_main(34),runoff
    
!   update
    real Npool(8),MineralN,NSN,NSNcrt1
    
!   output
    real NdemandT,CN(8),NresorpT
    real Nuptake,Nfix
    real NimmobT
    real Nscal_vcmax,Nscal_autoresp,Ccost_Nacq
    
!   internal
    real Ndeficit,outN(8),frac_soc(10),Nresorp(4),Nresorp_prop(4)
    real Creuse0,fwood_prop,Ccost_Nresorp
    real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M
    real Qroot0,Nup0,CcostperNfix,CcostperNuptake,MinNcrt1,MinNcrt2
    real Ndenitr_ref,MineralNref,Nvol_scal
    real Ndemand(3),Nimmob(8)
    real Ccost_Nupt,Ccost_Nfix,Nvol
    integer i,j
    real Nleach,Nminer,Nrunoff,water_leach
    
    
!    Nitrogen part ------------------------------------------------------------------
!   N deposit --------------------------------------------------------
    MineralN=MineralN+Ndeposit
  
!   Demand for producing new tissue (leaf, wood, and root) -------
    Ndemand(1)=NPP(1)/CN(1)+Cpool(1)/CNini(1)-Cpool(1)/CN(1)
    Ndemand(2)=NPP(2)/CN(2)
    Ndemand(3)=NPP(3)/CN(3)+Cpool(3)/CNini(3)-Cpool(3)/CN(3)
    NdemandT=Amax1(Ndemand(1)+Ndemand(2)+Ndemand(3),0.0)
    Ndeficit=Ndeficit+NdemandT

    Npool(1)=Npool(1)+Ndemand(1)
    Npool(2)=Npool(2)+Ndemand(2)
    Npool(3)=Npool(3)+Ndemand(3)

!    print*,"Ndeficit:",Ndeficit,NdemandT
    
!   N turnover out ---------------------------------------------------
    do i=1,8
        if(CN(i) .ne. 0)then  !edited by Chris
            OutN(i) = OutC(i)/CN(i)
            Npool(i)=Npool(i)-OutN(i)
        else
            OutN(i) = 0.0
        end if
    enddo
    
!   Part of plant outN is resorbed ---------------------------------
    Nresorp_prop(1)=0.6   ! Proportion of N resorbed before the dead of plant part
    Nresorp_prop(2)=0.5
    Nresorp_prop(3)=0.5
    Nresorp_prop(4)=0.6
    Creuse0=2.     ! C cost per N for resorption
    fwood_prop=0.15        ! 15% of woody litter is fine
    NresorpT=0.0
    Nresorp=0.0
    Ccost_Nresorp=0.0
    
    Nresorp(1)=OutN(1)*Nresorp_prop(1)*CN(1)/CNini(1)*NSNcrt1/NSN
!    print*,"Nresorp(1):",Nresorp(1),OutN(1),Nresorp_prop(1),CN(1),CNini(1),NSNcrt1,NSN
    Nresorp(2)=fwood_prop*OutN(2)*Nresorp_prop(2)*CN(2)/CNini(2)*NSNcrt1/NSN
    Nresorp(3)=(1-fwood_prop)*OutN(2)*Nresorp_prop(3)*CN(2)/CNini(2)*NSNcrt1/NSN
    Nresorp(4)=OutN(3)*Nresorp_prop(4)*CN(3)/CNini(3)
    NresorpT=Nresorp(1)+Nresorp(2)+Nresorp(3)+Nresorp(4)
    Ccost_Nresorp= Creuse0*NresorpT
    
!    print*,"Ccost_Nresorp:",Ccost_Nresorp, Creuse0,NresorpT,Nresorp(1),Nresorp(2),Nresorp(3),Nresorp(4)
    
    if(NresorpT.lt.Ndeficit)then
        Ndeficit=Ndeficit-NresorpT
    else
        NSN=NSN+(NresorpT-Ndeficit)
        Ndeficit=0.0
    endif
    
!   The other part of plant out N is transfered to litter -------------
    Npool(4)=Npool(4)+OutN(1)-Nresorp(1)+fwood_prop*OutN(2)-Nresorp(2)+OutN(3)-Nresorp(4)
    Npool(5)=Npool(5)+(1-fwood_prop)*OutN(2)-Nresorp(3)
    
!   Part of litter N and soil N are recycled through mineralization -------
    f_F2M=par_main(26)  ! initial: 0.55   
    f_C2M=par_main(27)  ! initial: 0.275    
    f_C2S=par_main(28)  ! initial: 0.275   
    f_M2S=par_main(29)  ! initial: 0.3
    f_M2P=par_main(30)  ! initial: 0.1
    f_S2P=par_main(31)  ! initial: 0.2
    f_S2M=par_main(32)  ! initial: 0.5
    f_P2M=par_main(33)  ! initial: 0.45
    
    Nminer=OutN(4)* (1. - f_F2M)   &
        &   +OutN(5)* (1. - f_C2M - f_C2S) &
        &   +OutN(6)* (1. - f_M2S - f_M2P) &
        &   +OutN(7)* (1. - f_S2P - f_S2M) &
        &   +OutN(8)* (1. - f_P2M)
        
    MineralN=MineralN+Nminer
    Npool(6)=Npool(6)+f_F2M*OutN(4)+f_C2M*OutN(5)+f_S2M*OutN(7)+f_P2M*OutN(8)
    Npool(7)=Npool(7)+f_C2S*OutN(5)+f_M2S*OutN(6)
    Npool(8)=Npool(8)+f_M2P*OutN(6)+f_S2P*OutN(7)
    
!   Further income of plant N through N uptake and fixation ----------------
!   New N is from uptake or biological fixation depend on the C cost of the processes
    Qroot0=100. ! unit: g/m2; initial: 500; may be adjusted
    Nup0 =0.02     ! nitrogen uptake rate
    CcostperNfix=12.      ! C cost per N for fixation
    CcostperNuptake=0.05      ! C cost per N for uptake
    MinNcrt1=0.10  ! Minimum mineral N content; gN/m2
    Nuptake=0.0
    Nfix=0.0
    Ccost_Nupt=0.0
    Ccost_Nfix=0.0
    
    If(Ndeficit>0.0)then
!        print*,"NSN 0:",NSN
        if(Ndeficit.lt.(NSN-NSNcrt1))then
            NSN=NSN-Ndeficit
!            print*,"NSN11:",NSN,Ndeficit
        else
            Ndeficit=Ndeficit-(NSN-NSNcrt1)
            NSN=NSNcrt1
!            print*,"NSN12:",NSN,Ndeficit,NSNcrt1
        endif
        
        if((MineralN.gt.MinNcrt1).and.(CcostperNuptake/MineralN<CcostperNfix))then  !  N uptake
!        print*,"test",MineralN,MinNcrt1,CcostperNuptake,CcostperNfix
            Nuptake=Amax1(AMIN1(Ndeficit,  &
                &   MineralN*Cpool(3)/(Cpool(3)+Qroot0), &
                &   Nup0/(CcostperNuptake/MineralN)),0.0)*NSNcrt1/NSN
            Ccost_Nupt=Nuptake*(CcostperNuptake/MineralN)
!            print*,"Ccost_Nupt:",Ndeficit,Nuptake,CcostperNuptake
            
            if(Nuptake.lt.Ndeficit) then
                MineralN=MineralN-Nuptake
                NSN=NSN
                Ndeficit=Ndeficit-Nuptake
            else
                MineralN=MineralN-Nuptake
                NSN=NSN+(Nuptake-Ndeficit)
                Ndeficit=0.0
            endif
            
!            print*,"Ndeficit",Ndeficit
            
        endif
!        print*,"NSN:",NSN,Ndeficit,MineralN,MinNcrt1
        if((NSN<24.*Ndeficit))then    ! Biological N fixation ! initial: * 30.
!        print*,"test",NSN,Ndeficit
            Nfix=Ndeficit
            Ccost_Nfix=CcostperNfix*Nfix
            if(Nfix.lt.Ndeficit)then
                Ndeficit=Ndeficit-Nfix
            else
                NSN=NSN+(Nfix-Ndeficit)
                Ndeficit=0.0
            endif
        endif
    endif
    
!    print*,"Ccost_Nfix:",Ccost_Nfix,CcostperNfix,Nfix,Ndeficit
!   print*,"Ccost_Nacq5:",Ccost_Nfix,CcostperNfix,Nfix,Ndeficit
   
!   Nitrogen immobilization -------------------------------------------
    Nimmob=0.
    NimmobT=0.
    
    if(MineralN>0)then
        do i=4,8
            if(CN(i).gt.CNini(i))  Nimmob(i)=Amax1(Cpool(i)/CNini(i)-Cpool(i)/CN(i),mineralN*0.1)
            Npool(i)=Npool(i)+Nimmob(i)
            if((MineralN-MinNcrt1).gt.Nimmob(i))then
                MineralN=MineralN-Nimmob(i)
            else
                MineralN=MinNcrt1   ! A minimum value
                Ndeficit=Ndeficit+Nimmob(i)-(MineralN-MinNcrt1)
            endif
        enddo
    endif
    
    NimmobT=Nimmob(4)+Nimmob(5)+Nimmob(6)+Nimmob(7)+Nimmob(8)
!    print*,"Ccost_Nacq4:",Ccost_Nfix
!   N runoff ------------------------------------------------
    MinNcrt2=1.00  ! Reference mineral N content; gN/m2
    
    if (MineralN.gt.MinNcrt1) then
        Nrunoff=runoff*1e-6*MineralN/MinNcrt2   ! unit of runoff: mm/h=kg/m2/h, assuming N conc. of 1 mg/L
        MineralN=MineralN-Nrunoff
    else
        Nrunoff=0.0
    endif
    
!   N leaching ---------------------------------------------
    if (MineralN.gt.MinNcrt1) then
        Nleach=water_leach*0.5e-6*MineralN/MinNcrt2
        MineralN=MineralN-Nleach 
    else
        Nleach=0.0
    endif
    
!   N volatilization ---------------------------------------
    frac_soc=(/0.1,0.2,0.25,0.15,0.12,0.08,0.04,0.02,0.02,0.02/)
    Ndenitr_ref=0.0002
    MineralNref=2.0 ! unit of g/m2
    Nvol_scal = 0.0 
!    print*,"Ccost_Nacq3:",Ccost_Nfix
    do j=1,10
        Nvol_scal = Nvol_scal+frac_soc(j)*Ndenitr_ref*Amax1(Amin1(exp((stemp(j)-25.)/10.),100.0),0.01)
    enddo
    
    Nvol=MineralN*Amin1(Nvol_scal*exp(MineralN/MineralNref-1),0.1)
    
    if (MineralN.gt.MinNcrt1) then
        MineralN=MineralN-Nvol 
    else
        Nvol=0.0
    endif
    
  
!   Update C/N ratio --------------------------------------
    do j=1,10
        if(Npool(j).eq.0)then
            CN(i)=CNini(i)
        else
            CN(i)=Cpool(i)/Npool(i)
        endif
    enddo
  
!    print*,"Npool:",Npool
!    print*,"Ccost_Nacq2:",Ccost_Nfix
!   N scalers for other subroutines ------------------------
!   If CN(i) increases, such as from 150 to 154, then Nscal_vcmax will be around 0.02
    Nscal_vcmax=exp(-CNini(1)*(CN(1)-CNini(1))/CNini(1)) ! /CNini(1) ! Duke
!    print*,"Nscal_vcmax:",Nscal_vcmax,CN(1),Cpool(1),Npool(1)
    
    Nscal_autoresp =exp(-(CN(1)-CNini(1))/CNini(1)) 

!   Total C cost for acquiring nitrogen -------------------
    Ccost_Nacq=0.0
    
    Ccost_Nacq=Ccost_Nupt+Ccost_Nfix+Ccost_Nresorp
    
!    print*,"Ccost_Nacq:",Ccost_Nacq,Ccost_Nupt,Ccost_Nfix
    
    if (isnan(MineralN)) then
        print*,"Warning: Mineral N is NAN !!!"
        stop
    endif
    
    return
end    
!   End of Eco_N subroutine =============================================================== 


!   ========================================================================================
!     carbon transfer according to Xu et al. 2007 
subroutine hetero_resp(vegetype,outC,tauC,Growth,litterfall,rootfall,radsol,par_main,Q10,   & ! input
                &   swater_slw,stemp,basew4sresp,   &   ! input
                &   Cpool,  &   ! update
                &   Rh_pools)   ! output
    
!   Define variables
    implicit none
!   input
    integer vegetype
    real outC(8),tauC(8),Growth(3),litterfall,rootfall,radsol,par_main(34),Q10
    real swater_slw,stemp(10),basew4sresp
    
!   update
    real Cpool(8)

!   output
    real Rh_pools(5)
    
!   internal
    real decom_wscal,decom_photo(2),decom_tscal(8),fwood_prop,Qplant
    real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M
    integer i
    character(len=80) commts

 
    Rh_pools=0.0
    
!   Plant carbon turnover -------------------------------------------------
    OutC(1)=litterfall
    OutC(2)=Cpool(2)/tauC(2)
    OutC(3)=rootfall ! Cpool(3)/tauC(3) ! *swater_slw
    
!    print*,"root:",rootfall
    
!   Heterotrophic respiration --------------------------------------------
    decom_photo(1)=0.0  ! 0.03/8760 ! according to Adair et al., 2017, Ecosphere.
    decom_photo(2)=0.01/8760 ! according to Adair et al., 2017, Ecosphere.
    
    decom_wscal=Amin1(swater_slw+basew4sresp,1.0)  ! Water scaler for decomposition ! May be adjusted later
    
    do i=4,8
        decom_tscal(i)=Amax1(Amin1(Q10**((stemp(2)-20.)/10.),1000.0),0.001)  ! Temperature scaler ! May be adjusted later
    enddo 
    
    OutC(4)=Cpool(4)/tauC(4)*decom_wscal* decom_tscal(4)+Cpool(4)*decom_photo(1)*radsol
    OutC(5)=Cpool(5)/tauC(5)*decom_wscal* decom_tscal(5) ! +Cpool(5)*decom_photo2*radsol
    OutC(6)=Cpool(6)/tauC(6)*decom_wscal* decom_tscal(6)
    OutC(7)=Cpool(7)/tauC(7)*decom_wscal* decom_tscal(7)
    OutC(8)=Cpool(8)/tauC(8)*decom_wscal* decom_tscal(8)

    
    
    f_F2M=par_main(26)  ! initial: 0.55   
    f_C2M=par_main(27)  ! initial: 0.275    
    f_C2S=par_main(28)  ! initial: 0.275   
    f_M2S=par_main(29)  ! initial: 0.3
    f_M2P=par_main(30)  ! initial: 0.1
    f_S2P=par_main(31)  ! initial: 0.2
    f_S2M=par_main(32)  ! initial: 0.5
    f_P2M=par_main(33)  ! initial: 0.45
    
!    print*,par_main,f_F2M

    Rh_pools(1)=OutC(4)* (1. - f_F2M)
    Rh_pools(2)=OutC(5)* (1. - f_C2M - f_C2S)
    Rh_pools(3)=OutC(6)* (1. - f_M2S - f_M2P)
    Rh_pools(4)=OutC(7)* (1. - f_S2P - f_S2M)
    Rh_pools(5)=OutC(8)* (1. - f_P2M)
    
!    print*,"Rh_pools",Rh_pools(1),outC(4),f_F2M,par_main(1)
!    read(*,*) commts
    
!   Update ecosystem C pools ------------------------------
    fwood_prop=0.15        ! 15% of woody litter is fine
    
 
    Cpool(1)=Amax1(Cpool(1)-OutC(1)+Growth(1),0.0)
!    print*,"Cpool(1):",Cpool(1),OutC(1),Growth(1)
    
    Cpool(2)=Amax1(Cpool(2)-OutC(2)+Growth(2),0.0)
    if(vegetype.eq.1) Cpool(2)=0.0
    
    Cpool(3)=Amax1(Cpool(3)-OutC(3)+Growth(3),0.0)
    Cpool(4)=Cpool(4)-OutC(4)+OutC(1)+fwood_prop*OutC(2)+OutC(3)
    Cpool(5)=Cpool(5)-OutC(5)+(1.-fwood_prop)*OutC(2)
    Cpool(6)=Cpool(6)-OutC(6)+f_F2M*OutC(4)+f_C2M*OutC(5)+ f_S2M*OutC(7)+f_P2M * OutC(8)
    Cpool(7)=Cpool(7)-OutC(7)+f_C2S*OutC(5)+f_M2S*OutC(6)
    Cpool(8)=Cpool(8)-OutC(8)+f_M2P*OutC(6)+f_S2P*OutC(7)
!    Q_plant =Cpool(1)+Cpool(2)+Cpool(3)
        
    return
end    
!   End of carbon_turnover subroutine =============================================================== 


!      ************************************************************************************************
!      *****************************   subroutines from methane and soil thermal modules **************
!      ************************************************************************************************

subroutine snow_d(lat,yr,days,  &  ! input: space and time
                &   rain_d,tair_dmean,  &   ! input: rain and temperature
                &   fa,fsub,rho_snow,decay_m,   &   ! input: snow related parameters
                &   dcount,snow_depth,melt)  ! outputs and/or inputs

    implicit none
!   inputs
    real lat
    integer yr,days
    real rain_d,tair_dmean
    real fa,fsub,decay_m
    real rho_snow ! snow density
!   outputs
    real dcount ! snow age in days
    real snow_depth,melt
!   internal variables
    real tr,dec,daylength  ! length the day time in the whole day
    real sublim,dsnow,snow_in
    real snow_dsim_pre
    real esat


!   Calculate daylength ------------------------------------------------------
    tr=0.0174532925
    dec=sin(((real(days)-70.)/365.)*360.*tr)*23.44
    daylength=acos(-tan(lat*tr)*tan(dec*tr))/7.5 
    daylength=daylength/tr/24.
    
!   Calculate snow age -------------------------------------------------------
    if (snow_depth .gt. 0.) then
        dcount = dcount +1.
    else 
        dcount =0.
    endif

!   Calculate sublimation rate and melting rate ------------------------------
    sublim=0.
    melt=0.
    if (tair_dmean .gt. 0. .and. snow_depth .gt. 0.) then
        sublim=fsub*715.5*daylength*(esat(tair_dmean)*0.001)/(tair_dmean+273.2)   ! 0.001 from Pa to kPa
        melt=fa*(2.63+2.55*tair_dmean+0.0912*tair_dmean*rain_d)   !dbmemo updated version
    endif
!   Consider snow age effect on melting --------------------------------------
    if (dcount .gt.0. .and. tair_dmean .lt.5.) then
        melt=melt*EXP(-decay_m*dcount/365.)  !dbmemo dbice
    endif
!   Assign precipitation as snow or rain -------------------------------------
    if (tair_dmean .le. 0.) then      
        snow_in =rain_d
    else
        snow_in = 0.
    endif
!   Update snow amount(dsnow) and snow depth(snow_dsimu; unit of mm) ----------
    dsnow=snow_in-sublim-melt  ! unit of mm/d
    snow_dsim_pre = snow_depth
    snow_depth =snow_depth+dsnow/rho_snow  ! unit of snow_depth: mm

    if (snow_depth .le. 0.0) then 
        snow_depth=0.0 
        melt = snow_dsim_pre*rho_snow +snow_in-sublim    !! for water part
    endif 
    melt=AMAX1(melt, 0.)  
    return
end ! subroutine snow_d
      
!    !    ========================================================================================
subroutine soil_temp(Rsoilab,RsoilabT,Esoil,&  ! input
                &   Tair,Cmolar,fbeam,LAI,sthick,layern,frlen,&  ! input
                &   wsmax,wsmin,swater_sld,swc,&  ! input
                &   parval_stem,Esoil_res,&  ! input
                &   sigma,emsoil,Rsoil, &  ! input
                &   Rconst,Rhocp,psyc,slope, & ! input: constants
                &   cpair,Patm,AirMa,H2OMw,H2OLv0,H2OLv,raero,zwt,&  ! input
                &   snow_depth,Tsnow, & ! input
                &   ice,waterv,ice_tw,&  ! both input and output
                &   wrt_tsoil,&   ! input
                &   Hsoil,G,stemp,sftmp,Troot) ! outputs 

    implicit none
    
!   Switches
    logical,parameter:: tsol_diagn=.false.
!   inputs
    integer layern
    real Rsoilab(3),RsoilabT,Esoil
    real Cmolar
    real Tair,fbeam,LAI,sthick(10),frlen(10)
    real wsmax(10),wsmin(10),swater_sld(10),swc(10)
    real parval_stem(20),sftmp,Esoil_res
    real Rconst,psyc
    real sigma,emsoil,cpair,Patm,AirMa,H2OMw,H2OLv0,H2OLv,raero,zwt
    real snow_depth,Tsnow,ice(10),waterv(10),ice_tw
    logical wrt_tsoil
    
!   outputs
    real Hsoil,G 
    real stemp(10),Troot
    
!   internal variables
    real tsoil,Twater
    real shcap_snow,condu_snow,condu_b,& 
    &   depth_ex,diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,b_bound

    real esat,theta_sat_min,Rsoil,difsv2,difsv1
    real delta
    real TairK
    real temph1,temph2     
    real hitmax,rflen,zopnd,thkns1,thkns2
    real flux_water,flux_snow
    real condu_water,shcap_water
    real albedo_water,heat_excess,heat_adjust
    real inter_var,latent_heat_fusion,QLsoil
    real rhocp,slope,fw1,resoil,rLAI,resdh,dnr,dsh,dgh,dle,drsdh
    real f_om,theta_sat_om,b_om,b_min,phi_om,phi_min,theta_sat,b_tot,phi_sat,gravi
    real water_table_depth,temph_water,temph_snow
    real condu_air,shcap_air,condu(10), shcap(10), tsoill_pre, thd_t
    real condu_soil,shcap_soil(10) 
    real resht_lai
    real condu_s,diff_air,d_cor,dcount
    real sftmp_pre
    
    real diff_water,water_tw
    real diff_ice,condu_ice,ice_density,ice_incr,shcap_ice,Tice
    
    integer n_layers,i

    real, allocatable ::depth_z(:) 
    allocate(depth_z(10))

    
!   Value assignments ---------------------------------------------------------------
    

    shcap_snow    =parval_stem(1)  ! tuneice worker better
    condu_snow    =parval_stem(2)
    condu_b       =parval_stem(3)  ! 
    depth_ex      =parval_stem(4)
    albedo_snow   =parval_stem(5)
    resht         =parval_stem(6) ! Initial: 40.

    shcap_soil    =parval_stem(11) ! EH mod; Unit: J/m3/k
 
    condu_soil    =parval_stem(12)  ! EH mod; unit: W/m/k
    shcap_water   =parval_stem(13) ! Unit: J/m3/k; Similar as in CLM5.0.
    condu_water   =parval_stem(14) ! Similar as in CLM5.0.
    shcap_air     =parval_stem(15)  ! Similar as in CLM5.0.
    condu_air     =parval_stem(16)   ! Similar as in CLM5.0.
    shcap_ice     =parval_stem(17)
    condu_ice     =parval_stem(18)
    latent_heat_fusion=parval_stem(19)
    ice_density   =parval_stem(20)
    
    drsdh=0.0
        
    thkns1=sthick(1)/2.
    thd_t=0.0

    diff_snow=3600.*condu_snow/shcap_snow*10000.
    diff_s=3600.*condu_b/shcap_soil(1)*10000.

    diff_water=3600.*condu_water/shcap_water*10000.
    diff_ice=3600.*condu_ice/shcap_ice*10000.
    
    diff_air=3600.*condu_air/shcap_air*10000. 
    
    water_tw=zwt*0.001-ice_tw ! might means total water that is liquid, add up all layers
    water_table_depth=zwt*0.1  
    
    
    if (water_table_depth .lt. 4. .and. water_table_depth .gt. 0.0) water_table_depth =0.    ! avoid numerical issues when 
    albedo_water =0.1      
  
    flux_snow = 0.0 
    depth_z=(/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 

   ! Total radiation absorbed by soil
    if (snow_depth .gt. 0.0) then
!        print*,"test1"
        RsoilabT=(Rsoilab(1)+Rsoilab(2))*(1-albedo_snow)/(1-0.1)+Rsoilab(3)  
    elseif (water_table_depth .gt. 0.0) then 
!        print*,"test2"
        RsoilabT=(Rsoilab(1)+Rsoilab(2))*(1-albedo_water)/(1-0.1)+Rsoilab(3)  
    else
        RsoilabT=Rsoilab(1)+Rsoilab(2)+Rsoilab(3)
    endif

!    print*,"RsoilabT:",RsoilabT,Rsoilab
    
    
!    thermodynamic parameters for air
!    H2OLv=H2oLv0-2.365e3*Tair
    slope=(esat(Tair+0.01)-esat(Tair))/0.01   

  
    fw1=AMIN1(AMAX1(1-swater_sld(1),0.3),1.0)
!      
    if (water_table_depth .gt. 0.0) Rsoil = 0. 

    rLAI=exp(LAI)     
!   aerodynamic resistance (s/m)
    resht_lai=Max1(resht*LAI,50.0)
!    print*,resht,resht_lai
    
!   Sensible heat
    Hsoil=rhocp*(sftmp-Tair)/resht_lai

    condu(1)=(1-wsmax(1))*condu_soil+(wsmax(1)-swc(1))*condu_air+swc(1)*condu_water
    shcap(1)=(1-wsmax(1))*shcap_soil(1)+(wsmax(1)-swc(1))*shcap_air+swc(1)*shcap_water
    difsv1=3600.*condu(1)/shcap(1)*10000.
    G=condu(1)*(sftmp-stemp(1))/(sthick(1)/2.*0.01)

!    print*,diff_s,difsv1
    
    if (snow_depth .gt. 0.0) then 
        G=condu_snow*(sftmp-Tsnow)/(snow_depth*0.1/2.)
    endif

!   Residual heat energy.
    RESDH=RsoilabT-Hsoil-Esoil-G
!   First derivative of net radiation; sensible heat; ground heat;
    DNR=4.*emsoil*sigma*(sftmp+273.2)**3   ! Correspond to RsoilabT
!    print*,"DNR:",DNR,emsoil,sigma,sftmp
    
    DSH=rhocp/resht_lai 	! Correspond to Hsoil
    
!    print*,"DSH:",DSH,rhocp,resht_lai 
    
    
    DGH=condu(1)/(sthick(1)/2.*0.01)		! Correspond to G
!    DLE=(DNR+DGH)*slope/(slope+psyc*(rsoil/raero+1.)) ! Same as for Esoil     ! Correspond to Esoil  
    DLE=(DNR-DGH)*slope/(slope+psyc*(rsoil/raero+1.)) ! Same as for Esoil     ! Correspond to Esoil  
!    print*,"DLE:",DLE,DNR,DGH,slope,psyc,rsoil,raero
    
    drsdh=-DNR-DSH-DLE-DGH     ! Initial
!    print*,"drsdh:",drsdh,DNR,DSH,DLE,DGH  
    
    ! Calculate increment DELTA.
!    DELTA=(resdh+Esoil_res)/drsdh
    DELTA=resdh/drsdh
    
!    print*,"delta:",sftmp,delta,resdh,drsdh,RsoilabT,Hsoil,Esoil,G
    
    sftmp_pre=sftmp
    sftmp=sftmp-DELTA
    
    if(sftmp.gt.60)then
        print*,"sftmp:",sftmp,delta,resdh,drsdh
        print*,"Warning: unreasonable surface temperature is calculated!!!"
        stop
    endif
    
    if (ABS(sftmp_pre -sftmp) .gt. 20. ) sftmp=sftmp_pre
    
    if(tsol_diagn) print*,"tsol_diag 3:",G,sftmp,DELTA,resdh,RsoilabT,Esoil
    
    do i=1,layern
        if(tsol_diagn) print*,"tsol_diag 4:",i,stemp(i)
        
        Tsoill_pre=stemp(i)
        
        if (water_table_depth .lt. 0.0 .and. -water_table_depth .lt. depth_z(i)) then
            waterv(i)=wsmax(i)-ice(i)
        else
            waterv(i)=swc(i)-ice(i)
        endif
    
        if (i .eq. 1) then 
            depth_z(1)=sthick(1)
        else 
            depth_z(i)=depth_z(i-1)+sthick(i)
        endif

        thkns2=(sthick(i)+sthick(i+1))/2.

        if (i .eq. layern) then
            difsv2=3600.*condu(i)/shcap(i)*10000. 
        else
            condu(i+1)=(1-wsmax(i+1))*condu_soil+(wsmax(i+1)-swc(i+1))*condu_air+ &
            &   waterv(i+1)*condu_water+ice(i+1)*condu_ice
            shcap(i+1)=(1-wsmax(i+1))*shcap_soil(i+1)+(wsmax(i+1)-swc(i+1))*shcap_air+ &
            &   waterv(i+1)*shcap_water+ice(i+1)*shcap_ice
            
            difsv2=3600.*condu(i+1)/shcap(i+1)*10000.
        endif
        
        temph2=(difsv1+difsv2)*(stemp(i)-stemp(i+1))/thkns2 

        
!!!!!!!!!!!!!!!!!!!! start first layer !!!!!!!!!!!!!!!!!!!!!!
        !!!!!! adjust if there are snow or water layer above !!!!!!!!!!!!!!!!!!!!
        if(i.eq.1) then
            if (snow_depth .gt. 0.) then   
                print*,"Snow depth is great than 0."
                temph_snow = Amin1(diff_snow,difsv1)*(Tsnow-stemp(1))/((snow_depth*0.1+sthick(1))/2.)
                Tsnow=Tsnow+(exp(-depth_ex*snow_depth*0.1)*diff_snow*(sftmp-Tsnow)/(snow_depth*0.1/2.) &
                    &   -temph_snow)/(snow_depth*0.1/2.+(snow_depth*0.1+sthick(1))/2.) 

                stemp(1)=stemp(1)+(temph_snow-temph2)/((snow_depth*0.1+sthick(1))/2.+thkns2) 

                if (Tsnow .gt.0.0) then 
                    Tsnow =0.0   
                    stemp(1)=0.
                endif
                
                drsdh =0.0    ! temporarily set drsdh =0 for heat adjustment of soil when  
                sftmp= (stemp(1)+Tsnow)/2.
                
            elseif (water_table_depth .gt. 0.) then  
!                print*,"Water table depth is greater than 0."
                temph_water=(diff_water+difsv1)*(Twater-stemp(1))/((water_table_depth+sthick(1))/2.)! there is snow layer 
                Twater=Twater+(2.*diff_water*(sftmp-Twater)/(water_table_depth/2.-temph_water))/ &
                    &   (water_table_depth/2.+(water_table_depth+sthick(1))/2.) 
                    
!               Phase change of surface water 
                if (Twater .lt. 0.0 .and. water_tw .gt. 0.0) then  ! freeze 
                    heat_excess=-(shcap_water/3600*water_tw-drsdh)*Twater
                    ice_incr=heat_excess*3600./latent_heat_fusion/ice_density
                    if (ice_incr .lt. water_tw) then
                        ice_tw=ice_tw +ice_incr
                        water_tw=water_tw-ice_incr
                        Twater=0.0
                        Tice=0.0
                    else
                        ice_tw=ice_tw +water_tw
                        water_tw=0.0
                        Tice = Tice - latent_heat_fusion*(ice_incr-water_tw)*ice_density/(shcap_ice*ice_tw)
                    endif     
                elseif (Twater .gt. 0.0 .and. ice_tw .gt. 0.0) then    ! thraw              
                    heat_excess=(shcap_water/3600*ice_tw-drsdh)*Twater
                    ice_incr=heat_excess*3600./latent_heat_fusion/ice_density
                    if (ice_incr .lt. ice_tw) then
                        ice_tw=ice_tw-ice_incr
                        water_tw=water_tw+ice_incr
                        Twater=0.0
                        Tice=0.0
                    else
                        water_tw=water_tw +ice_tw
                        ice_tw=0.0
                        Twater = Twater + latent_heat_fusion*(ice_incr-ice_tw)*ice_density/(shcap_water*water_tw)
                    endif
                endif         
                    
                temph2=(difsv1+diff_water)*(stemp(i)-stemp(i+1))/thkns2

                if (water_tw .eq. 0.0 .and. ice_tw .gt. 0.0) then 
                    stemp(1)=stemp(1)+(2.*diff_ice*(Tice-stemp(1))/thkns1-temph2)/(thkns1+thkns2) 
                else 
                    stemp(1)=stemp(1)+(2.*diff_water*(Twater-stemp(1))/thkns1-temph2)/(thkns1+thkns2) 
                endif    

                drsdh =0.0    ! temporarily set drsdh =0 for heat adjustment of soil  
                
            else   ! without snow or high (than surface) water table 
                stemp(1)=stemp(1)+(diff_s*(sftmp-stemp(1))/thkns1-temph2)/(thkns1+thkns2)
            endif
!            print*,"stemp(1) 1:",stemp(1)
            
!           Phase change in top soil ------------------------------------------------------------ 
            heat_excess=drsdh*(thd_t-stemp(i))+shcap(i)*sthick(i)*(stemp(i)-thd_t)/360000.
            ice_incr=heat_excess*3600./latent_heat_fusion/ice_density
            inter_var = ice(i)   
            if (ice_incr .lt. 0.) then     ! freeze             
                ice(i)=Amin1(waterv(i)+inter_var,ice(i)-ice_incr/(sthick(i)*0.01))            
            else 
                ice(i)= Amax1(ice(i)-ice_incr/(sthick(i)*0.01),0.0)              
            endif
            heat_adjust=heat_excess-latent_heat_fusion*((inter_var-ice(i))*sthick(i)*0.01)*ice_density/3600.
            stemp(i)=thd_t+heat_adjust/(shcap(i)*sthick(i)/360000.-drsdh) 
            
!            print*,"stemp(1) 2:",i,heat_adjust,heat_excess,drsdh,thd_t,stemp(i),shcap(i),sthick(i),stemp(i),thd_t
            
!            if(tsol_diagn) print*,"tsol_diag 5:",diff_s,sftmp,stemp(1)
            
        else  ! after this for soil layers 2-10
            if ( i .gt. (layern-1)) then 
                temph2=0.00003
                thkns2=500  ! boundary conditions, rethink
            endif            

            stemp(i)=stemp(i)+(temph1-temph2)/(thkns1+thkns2) 
            heat_excess=shcap(i)*sthick(i)*(stemp(i)-thd_t)/360000. 
            ice_incr=heat_excess*3600./latent_heat_fusion/ice_density                     
            inter_var = ice(i) 
            if (ice_incr .lt. 0.) then     ! freeze             
               ice(i)=Amin1(waterv(i)+inter_var,ice(i)-ice_incr/(sthick(i)*0.01))             
            else 
               ice(i) = Amax1(ice(i)-ice_incr/(sthick(i)*0.01),0.0)              
            endif         
            
            heat_adjust=heat_excess-latent_heat_fusion*((inter_var-ice(i))*sthick(i)*0.01)*ice_density/3600.
            stemp(i)=thd_t+heat_adjust/(shcap(i)/360000.*sthick(i))
            
!            print*,"tsol_diag 5_1:",stemp
            
        endif
        if (ABS(tsoill_pre -stemp(i)) .gt. 5. ) stemp(i)=tsoill_pre
!        print*,"tsol_diag 5_2:",i,stemp(i)
        
        if (isnan(stemp(i))) then
            print*,"soil temperature:",i,heat_excess,drsdh,thd_t,stemp(i)
            print*,"Warning: soil temperature is NAN!!!"
            stop
        endif
        
        if (stemp(i).gt.60 .or. stemp(i).lt.-40) then
            print*,"soil temperature:",i,stemp(i)
            print*,"Warning: extreme soil temperature is calculated!!!"
            stop
        endif
        
        TEMPH1=TEMPH2
        THKNS1=THKNS2
        DIFSV1=DIFSV2        
    enddo
    
    Troot=0.0
    do i=1,layern       
        Troot=Troot+stemp(i)*frlen(i)  ! linear (l) soil water scaler (s) weighted by the whole (w) root system 
    enddo

    deallocate(depth_z)
    
    if(tsol_diagn) print*,"tsol_diag 6:",stemp
    return 
end ! subroutine soil_temp
!   ***********************************************************************************

!   ============================================================================================
!   Plant respiration subroutine
subroutine plant_respiration(vegetype,Tleaf,sftmp,Troot,Cpool,GPP,LAI,swater_slw,gcostpro,mresp20,   &   ! input
                    &   Q10paccrate, & ! input
                    &   RmainT,Rmain,RgrowthT,Rgrowth)  ! output

    implicit none
!   input
    integer vegetype
    real Tleaf,sftmp,Troot,Cpool(8),GPP,LAI,swater_slw,gcostpro,mresp20(3),Q10paccrate(3)
!   output
    real RmainT,Rmain(3),RgrowthT,Rgrowth(3)
!   internal variables
    real Q10plant(3),LAI4resp
    
!   Plant respiration -------------------------------------------------------
    RmainT=0.0
    RgrowthT=0.0
    Rmain=0.0
    Rgrowth=0.0
    LAI4resp=1.0   
    
!   Plant maintenance respiration, according to Arora, V. K. (2003) 
!       Varied Q10, reflect the acclimation of Q10 to temperature
!       Set according to Atkin and Tjoelker, 2003        
    Q10plant(1)=Amax1(3.09-Q10paccrate(1)*sftmp,1.0)  ! may change to Tleaf later
    Q10plant(2)=Amax1(3.00-Q10paccrate(2)*sftmp,1.0)
    Q10plant(3)=Amax1(3.00-Q10paccrate(3)*Troot,1.0)

    Rmain(1)=mresp20(1)*Cpool(1)*Q10plant(1)**((sftmp-20)/10)*LAI4resp !*(1-swater_slw)
    Rmain(2)=mresp20(2)*Cpool(2)*Q10plant(2)**((sftmp-20)/10)*LAI4resp !*(1-swater_slw)
    if(vegetype.eq.1) Rmain(2)=0.0
!    print*,"Rmain(2):",vegetype,Rmain(2)
    
    Rmain(3)=mresp20(3)*Cpool(3)*Q10plant(3)**((Troot-20)/10)*LAI4resp !*(1-swater_slw)

    RmainT=sum(Rmain)

!    print*,"RmainT:",RmainT,Rmain,Cpool
    
!   Plant growth/construction respiration, according to Ryan, M. G. (1991)
    RgrowthT=Amax1(GPP-RmainT,0.0)*gcostpro*LAI4resp   ! gcostprop initial: 0.1
    
    return
end     ! End of plant_respiration subroutine

!     ========================================================================================
!     subroutines used by canopy submodel
!    the multi-layered canopy model developed by 
!    Ray Leuning with the new radiative transfer scheme   
!    implemented by Y.P. Wang (from Sellers 1986)
!    12/Sept/96 (YPW) correction for mean surface temperature of sunlit
!    and shaded leaves
!    Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)} 


subroutine xlayers(m,yr,yr_spinup,phenoset,  & ! input
            &  radabv,fbeam,co2ca,Dair,eairP,Tair,TairK,wind,sftmp,Vcmax_par,& ! input     ! 
            &  wsmax,wsmin,swater_slw,swater_sld,LAI, & ! input
            &  coszen,alpha,Rconst,rhocp,Cmolar,psyc,slope,& ! input
            &  tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,& ! input
            &  sigma,emleaf,emsoil,Ds0,& ! input
            &  cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat,& ! input
            &  gsc0,stom_n,raero,Rsoil,& ! input
            &  Vcmax_n,Jmax_n,Kc_ref,conKo0,Ekc,Eko,o2ci,& ! input
            &  gam0,gam1,gam2,& ! input
            &  G, & ! input
            &  RsoilabT,Rsoilab,&  ! output
            &  Acanop,Etransp,QLleaf,QLair,extKb,&  ! output
            &  Esoil,Eabsorp_leaf,Hcan,gsc_l,Tleaf)  ! output
    implicit none
   
!   Switches
    logical,parameter:: lay_diagn=.false.
    logical,parameter:: wrt_acanop=.false.
    
    real,external::esat
    
!   input variables
    integer m,yr,yr_spinup,phenoset
    real radabv(2),fbeam,co2ca,Dair,eairP,Tair,wind,sftmp,Vcmax_par(3)
    real wsmax(10),wsmin(10),swater_slw,swater_sld(10),LAI
    real coszen,alpha,slope
    real tauL(3),rhoL(3),rhoS(3),xfang,extkd,extkU,wleaf
    real sigma,emleaf,emsoil,Ds0
    real cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat
    real gsc0,stom_n
    real Vcmax_n,Jmax_n,Kc_ref,conKo0,Ekc,Eko,o2ci
    real gam0,gam1,gam2
    real G
    
!   output variables
    real RsoilabT   ! need check
    real Acanop,Etransp
    real Rsoilab(3)
    real QLleaf,QLair
    real Esoil,Tleaf
          
!   internal variable
    real Rconst
    real Acan(2),Aleaf(2)
    real Cmolar,cozen15,cozen45,cozen75,co2ci(2)
    real Etransp_side(2),Etransp_ng(2),emair,extkb,extkn
    real flait1,flai,funG,fw1,fshd,fslt
    real Gaussx(5),Gaussw(5),gbleaf(2),Gbwc1,Gbwc2,gsleaf(2),grn,Gswc1,Gswc2,gsc_l1,gsc_l2,gsc_l
    real Hcan1,Hcan2,Hleaf(2),Hcan
    real Jmax_nl,kpr(3,2)
    real pi180,psyc
    real Qabs(3,2), Qcan1,Qcan2,QLsoil
    real Rcan1,Rcan2,raero,Eabsorp_leaf,Eabsorp_leaf_side(2),Eabsorp_leaf_ng(2)
    real rhoc(3,2),reff(3,2),rhoc15,rhoc45,rhoc75,rhoch,rhocp,Rsoil
    real scatt(3)
    real TairK,Tleaf_ng(2),Tleaf_side(2),transd
    real Vcmax_nl,wind_l
    real xK15,xK45,xK75,xphi1,xphi2
    integer ng,nw
     
! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!* 5-point
    data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
    data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/

!     reset the variables
    xfang=0.0
    extkU=0.51
    
    
    Eabsorp_leaf_side=0.0  ! 1 for sunlit leaf; 2 for shaded leaf
    Qcan1=0.0        !vis rad
    Qcan2=0.0
    Rcan1=0.0        !NIR rad
    Rcan2=0.0
    Acan=0.0        !CO2
    Acanop=0.0
    Etransp_side=0.0
    Etransp=0.0
    Hcan1=0.0        !Sens heat
    Hcan2=0.0
    Gbwc1=0.0        !Boundary layer conductance
    Gbwc2=0.0
    Gswc1=0.0        !Canopy conductance
    Gswc2=0.0
    gsc_l1=0.0
    gsc_l2=0.0
    Tleaf_side=0.0

    Esoil=0
    
    co2ci=co2ca*0.9 ! set initial co2ci values

    call light_extinction_rate(xfang,coszen,Rhos,tauL,rhoL,LAI, &  ! input
                &   extKb,extkd,extkn,kpr,scatt,reff)  ! output  
    
                
    TairK=Tair+273.2

    if(lay_diagn) print*,"lay_diag 2:"
    
    do ng=1,5
        flai=gaussx(ng)*LAI
!        print*,"FLAI:",flai,gaussx(ng),LAI
!        radiation absorption for visible and near infra-red
        call Qabs_cal(flai,coszen,radabv,fbeam,reff,kpr,scatt,xfang,Qabs)
        
!        isothermal net radiation & radiation conductance at canopy top
        call Radiso(flai,LAI,Qabs,extkd,Tair,TairK,eairP,cpair,Patm,    & ! input
            &   fbeam,airMa,Rconst,sigma,emleaf,emsoil,rhocp,   &    ! input
            &   emair,Eabsorp_leaf_ng,grn)    ! output

        wind_l=wind*exp(-extkU*flai) ! wind speed adjusted by LAI
        Vcmax_nl=Vcmax_n*exp(-extkn*flai)
!        Amin1(swater_sld(1)+0.1,1.0) ! adjust Vcmax_n by LAI
!        print*,"Vcmax_nl:",Vcmax_nl,Vcmax_n
        
        Jmax_nl=Jmax_n*exp(-extkn*flai) ! adjust Jmax_n by LAI

        if(lay_diagn) print*,"lay_diag 3:"
!        print*,"Vcmax_nl:",Qabs,flai,coszen,radabv,fbeam,reff,kpr,scatt
        call one_layer(m,ng,yr,yr_spinup,phenoset,radabv,Qabs,Eabsorp_leaf_ng,grn,wind_l,Tair,TairK,Dair,      &    ! input
        &   co2ca,wleaf,raero,Ds0,swater_slw,  &    ! input
        &   Rconst,rhocp,cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat,  &    ! input
        &   gsc0,alpha,stom_n,Cmolar,psyc,slope,       &    ! input
        &   Vcmax_nl,Jmax_nl,Vcmax_par,          &    ! input
        &   gam0,gam1,gam2,         &    ! input
        &   Aleaf,Etransp_ng,Hleaf,Tleaf_ng,gbleaf,gsleaf,co2ci)    ! output

        if(lay_diagn) print*,"lay_diag 4:"
                        
        fslt=exp(-extKb*flai)                        !fraction of sunlit leaves
        fshd=1.0-fslt                                !fraction of shaded leaves
        Eabsorp_leaf_side(1)=Eabsorp_leaf_side(1)+fslt*Eabsorp_leaf_ng(1)*Gaussw(ng)*LAI  !Isothermal net radiation
        Eabsorp_leaf_side(2)=Eabsorp_leaf_side(2)+fshd*Eabsorp_leaf_ng(2)*Gaussw(ng)*LAI

        Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*LAI  !visible
        Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*LAI

        Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*LAI  !NIR
        Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*LAI

        if(Aleaf(1).lt.0.0)Aleaf(1)=0.0      !Weng 2/16/2006
        if(Aleaf(2).lt.0.0)Aleaf(2)=0.0      !Weng 2/16/2006
!       C gain
        Acan(1)=Acan(1)+fslt*Aleaf(1)*Gaussw(ng)*stom_n*LAI    !amphi/hypostomatous
        Acan(2)=Acan(2)+fshd*Aleaf(2)*Gaussw(ng)*stom_n*LAI   ! CO2
!       Energy for evaporation
        Etransp_side(1)=Etransp_side(1)+fslt*Etransp_ng(1)*Gaussw(ng)*LAI     ! Transpiration
        Etransp_side(2)=Etransp_side(2)+fshd*Etransp_ng(2)*Gaussw(ng)*LAI

        Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*LAI     ! Sensible heat
        Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*LAI

        Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*LAI*stom_n
        Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*LAI*stom_n

        Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*LAI*stom_n   ! Canopy conductance / Scale up
        Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*LAI*stom_n
        
        Gsc_l1=Gsc_l1+fslt*gsleaf(1)*Gaussw(ng)  ! leaf scale of co2; unit of mol/m2/s
        Gsc_l2=Gsc_l2+fshd*gsleaf(2)*Gaussw(ng)
        
!       Leaf temperature (oC)
        Tleaf_side(1)=Tleaf_side(1)+fslt*Tleaf_ng(1)*Gaussw(ng)*LAI       ! Leaf temperature
        Tleaf_side(2)=Tleaf_side(2)+fshd*Tleaf_ng(2)*Gaussw(ng)*LAI
        
!        print*,"Tleaf_side(2):",ng,Tleaf_side(2),fshd,fslt,Tleaf_ng(2),Gaussw(ng),LAI

    enddo  ! 5 layers
    
    if(lay_diagn) print*,"lay_diag 5:"
        
    Acanop=sum(Acan)  ! EH mod
!    print*,"Acanop:",Acanop
    
    Etransp=Etransp_side(1)+Etransp_side(2)  ! EH mod
    Hcan=Hcan1+Hcan2
    Eabsorp_leaf=Eabsorp_leaf_side(1)+Eabsorp_leaf_side(2)
    gsc_l=Gsc_l1+Gsc_l2
    
    FLAIT1=(1.0-exp(-extKb*LAI))/extkb
    Tleaf_side(1)=Tleaf_side(1)/FLAIT1
    Tleaf_side(2)=Tleaf_side(2)/(LAI-FLAIT1)
    
    Tleaf=(Tleaf_side(1)+Tleaf_side(2))/5   ! 5 layers
    
!    print*,"Tleaf:",Tleaf
    
    if(Tleaf.gt.60 .or. Tleaf.lt.-30) then
        print*,"Leaf temperature:",Tleaf,Tleaf_side(1),Tleaf_side(2),Tleaf_ng(2)
        print*,"Warning: unreseanable leaf temperature is calculated!!!"
!        stop
    endif
    
!   Soil surface energy and water fluxes
!   Radiation absorbed by soil
    Rsoilab(1)=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*LAI)        &
            &       +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*LAI)          !visible sorption coefficient
    Rsoilab(2)=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*LAI)        &
            &       +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*LAI)          !NIR sorption coefficient

    Rsoilab(1)=Rsoilab(1)*Radabv(1)  ! Radabv (1) equal to half of solar radiation  
    Rsoilab(2)=Rsoilab(2)*Radabv(2)  ! Radabv (2) equal to half of solar radiation  
 

    QLair=emair*sigma*(TairK**4)*exp(-extkd*LAI)
    QLleaf=emleaf*sigma*((Tleaf_side(1)+273.2)**4)*exp(-extkb*LAI)           &
        &   +emleaf*sigma*((Tleaf_side(2)+273.2)**4)*(1.0-exp(-extkb*LAI))
    QLleaf=QLleaf*(1.0-exp(-extkd*LAI))
    QLsoil=emsoil*sigma*((sftmp+273.2)**4)  ! longwave radiation out (W/m2)

    Rsoilab(3)=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil

!    print*,"Rsoilab3:",Rsoilab(3),(QLair+QLleaf),(1.0-rhoS(3)),QLsoil,emsoil,sigma,sftmp
   
!    Total radiation absorbed by soil    
    RsoilabT=SUM(Rsoilab)
!    print*,"RsoilabT:",RsoilabT,Rsoilab
    
    if (RsoilabT.gt.1400) then
        print*,"RsoilabT:",RsoilabT,QLleaf,Tleaf_side
        print*,"Warning: unreasonable high RsoilabT!!!"
        stop
    endif
    
    if(lay_diagn) print*,"lay_diag 6:"
    
!    thermodynamic parameters for air

!    H2OLv=H2oLv0-2.365e3*Tair
    slope=(esat(Tair+0.1)-esat(Tair))/0.1
    

!    fw1=AMIN1(AMAX1((wsmax(1)-swc(1))/(wsmax(1)-wsmin(1)),0.05),1.0)
  

!   According to Penman-Monteith equation, a reference: http://www.fao.org/docrep/X0490E/x0490e06.htm
!   RsoilabT: soil sorption of radiation/net radiation; 
!   G: ground heat
!   slope: unit of KPa/oC, represents the slope of saturation vapour pressure temperature relationship
!   Dair:  saturation vapour pressure deficit (kPa)
!   raero: aerodynamic resistance (s/m)
!   rsoil: bulk surface resistance (s/m)     
!   psyc: pyschometric constant; unit: kPa/oC
!   Esoil/RsoilabT/G have the same units
    Esoil=(slope*(RsoilabT-G)+rhocp*Dair/raero)/       &
        &  (slope+psyc*(rsoil/raero+1.))   ! Unit of KPa/oC
!    print*,"Esoil:",Esoil,RsoilabT,G,slope*(RsoilabT-G),RsoilabT,Tair,Tleaf           !!!!!!!!!!!!!!!!!!!!! YZhou    
!   Esoil=0.3*(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &  !0.3 added by Chris will adjusted later
!     &      (slope+psyc*(Rsoil/(raero+rLAI)+1.))
        
     if(yr.gt.yr_spinup .and. wrt_acanop) then
         write(3201,32011) m,swater_slw,Acanop,Acan(1),Acan(2),LAI
     endif
32011   format(1(I8,","),4(f15.10,","),(f15.10))
        
        
    return
end ! subroutine xlayers
!   =======================================================================================


!   Calculating light extinction rate ===================================================    
subroutine light_extinction_rate(xfang,coszen,Rhos,tauL,rhoL,LAI, &  ! input
                &   extKb,extkd,extkn,kpr,scatt,reff)  ! output  
                !    use the radiation scheme developed by
!    Goudriaan (1977, Goudriaan and van Larr 1995)
!=================================================================
!    Variable      unit      defintion
!    LAI         m2/m2     canopy leaf area index       
!    coszen                  cosine of the zenith angle of the sun
!    radabv(nW)    W/m2      incoming radiation above the canopy
!    fbeam                   beam fraction
!    fdiff                   diffuse fraction
!    funG(=0.5)              Ross's G function
!    extkb                   extinction coefficient for beam PAR
!    extkd                   extinction coefficient for diffuse PAR
!    albedo                  single scattering albedo
!    scatB                   upscattering parameter for beam
!    scatD                   upscattering parameter for diffuse
! ==================================================================
!    all intermediate variables in the calculation correspond
!    to the variables in the Appendix of of Seller (1985) with
!    a prefix of "x".
                
implicit none
!   input    
    real xfang,coszen,Rhos(3),tauL(3),rhoL(3),LAI
!   output
    real extKb,extkd,extkn
    real kpr(3,2),reff(3,2)
!   internal
    real xphi1,xphi2,funG
    real pi180,cozen15,cozen45,cozen75,xK15,xK45,xK75,transd
    real scatt(3),rhoc(2,2),rhoch,rhoc15,rhoc45,rhoc75
    integer nw
    

!    Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
    xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
    xphi2 = 0.877 * (1.0 - 2.0*xphi1)
    funG=xphi1+xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam

    if(coszen.gt.0) then    ! Day
        extKb=funG/coszen           !beam extinction coeff - black leaves
    else    ! Night
        extKb=100.             
    end if

!     Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
    pi180=3.1416/180.
    cozen15=cos(pi180*15)
    cozen45=cos(pi180*45)
    cozen75=cos(pi180*75)
    xK15=xphi1/cozen15+xphi2
    xK45=xphi1/cozen45+xphi2
    xK75=xphi1/cozen75+xphi2

    transd=0.308*exp(-xK15*LAI)+0.514*exp(-xK45*LAI)+0.178*exp(-xK75*LAI)
    extkd=(-1./LAI)*alog(transd)
    extkn=extkd                        !N distribution coeff 
    
!canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
!                                               second; 1=beam, 2=diffuse
    do nw=1,2                                      !nw:1=VIS, 2=NIR
        scatt(nw)=tauL(nw)+rhoL(nw)                      !scattering coeff
        if((1.-scatt(nw))<0.0)scatt(nw)=0.9999           ! Weng 10/31/2008
        kpr(nw,1)=extKb*sqrt(1.-scatt(nw))               !modified k beam scattered (6.20)
        kpr(nw,2)=extkd*sqrt(1.-scatt(nw))             !modified k diffuse (6.20)
        rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))       !canopy reflection black horizontal leaves (6.19)
        rhoc15=2.*xK15*rhoch/(xK15+extkd)            !canopy reflection (6.21) diffuse
        rhoc45=2.*xK45*rhoch/(xK45+extkd)
        rhoc75=2.*xK75*rhoch/(xK75+extkd)
        rhoc(nw,2)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
        rhoc(nw,1)=2.*extKb/(extKb+extkd)*rhoch          !canopy reflection (6.21) beam 
        reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))   &   !effective canopy-soil reflection coeff - beam (6.27)
        &            *exp(-2.*kpr(nw,1)*LAI) 
        reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))   &   !effective canopy-soil reflection coeff - diffuse (6.27)
        &            *exp(-2.*kpr(nw,2)*LAI)  
    enddo
    
    return
end ! subroutine light_extinction_rate
!   ==============================================================


!   =======================================================================================
subroutine Qabs_cal(FLAI,coszen,radabv,fbeam,reff,kpr,scatt,xfang, &    ! input
                    Qabs)                                     ! output

!    for spheric leaf angle distribution only
!    compute within canopy radiation (PAR and near infra-red bands)
!    using two-stream approximation (Goudriaan & vanLaar 1994)
!    tauL: leaf transmittance
!    rhoL: leaf reflectance
!    rhoS: soil reflectance
!    sfang XiL function of Ross (1975) - allows for departure from spherical LAD
!         (-1 vertical, +1 horizontal leaves, 0 spherical)
!    FLAI: canopy leaf area index
!    funG: Ross' G function
!    scatB: upscatter parameter for direct beam
!    scatD: upscatter parameter for diffuse
!    albedo: single scattering albedo
!    output:
!    Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
!                       type=1 for sunlit;   =2 for shaded (W/m2)
                    
    implicit none        
!   input variables
    real FLAI,coszen,radabv(2),fbeam,reff(3,2),kpr(3,2),scatt(2),xfang
    
!   output variables
    real Qabs(3,2)

!   internal variables
    real extKb,funG,Qb0,Qd0,xphi1,xphi2
    integer nw
    
!     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
    xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
    xphi2 = 0.877 * (1.0 - 2.0*xphi1)
    funG=xphi1+xphi2*coszen                            !G-function: Projection of unit leaf area in direction of beam

    if(coszen.gt.0) then                                  !check if day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
    else
        extKb=100.
    end if

! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
    do nw=1,2
        Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
        Qb0=fbeam*radabv(nw)                                               !beam incident radiation
        Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
            &   Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
            &   extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
        Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves
        
!        print*,"Qabs_internal:",nw,Qabs(nw,2),extKb,Qb0,scatt(nw)
        
        
    end do
    return
end ! subroutine goudriaan
!   *****************************************************************************
      
!****************************************************************************
subroutine Radiso(flai,LAI,Qabs,extkd,Tair,TairK,eairP,cpair,Patm,    & ! input
            &   fbeam,airMa,Rconst,sigma,emleaf,emsoil,rhocp,   &    ! input
            &   emair,Eabsorp_leaf_ng,grn)    ! output
!     output
!     Eabsorp_leaf_ng(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
!     23 Dec 1994
!     calculates isothermal net radiation for sunlit and shaded leaves under clear skies
!     implicit real (a-z)

    implicit none
!    input variables
    real flai,LAI,Qabs(3,2),extkd,Tair,TairK,eairP,cpair,Patm
    real fbeam,airMa,Rconst,sigma,emleaf,emsoil
    real rhocp
    
!    output variables
    real Eabsorp_leaf_ng(2)
    real emair,grn
    
!    internal variables
    real Bn0,Bnxi,emcloud,emsky,ep8z,tau8
   

! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
!   The coefficient was originally set as 0.642; however, it should be 
    emsky=0.642*(eairP/TairK)**(1./7)       !note eair in Pa 

! apparent emissivity from clouds (Kimball et al 1982)
    ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)
    tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
    emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10

! apparent emissivity from sky plus clouds      
    emair=emsky+emcloud

    if(emair.gt.1.0) emair=1.0

! net isothermal outgoing longwave radiation per unit leaf area at canopy
! top & thin layer at flai (Note Rn* = Sn+Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
    Bn0=sigma*(TairK**4.)
    Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
        & +exp(-extkd*(LAI-flai))*(emsoil-emleaf))
!     isothermal net radiation per unit leaf area for thin layer of sunlit and
!     shaded leaves
    Eabsorp_leaf_ng(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
    Eabsorp_leaf_ng(2)=Qabs(1,2)+Qabs(2,2)+Bnxi
    
!    print*,"Eabsorp_leaf_ng:",Qabs(1,1),Qabs(2,1)

!     radiation conductance (m/s)
!     According to D7 in Leuning et al. 1995
    grn=4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
        &    (exp(-extkd*flai)+exp(-extkd*(LAI-flai)))       &
        &    /rhocp
    return
end   ! subroutine Radiso
!   ******************************************************************************

!     ****************************************************************************
subroutine one_layer(m,ng,yr,yr_spinup,phenoset,radabv,Qabs,Eabsorp_leaf_ng,grn,wind_l,Tair,TairK,Dair,      &    ! input
        &   co2ca,wleaf,raero,Ds0,swater_slw,  &    ! input
        &   Rconst,rhocp,cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat,  &    ! input
        &   gsc0,alpha,stom_n,Cmolar,psyc,slope,               &    ! input
        &   Vcmax_nl,Jmax_nl,Vcmax_par,          &    ! input
        &   gam0,gam1,gam2,         &    ! input
        &   Aleaf,Etransp_ng,Hleaf,Tleaf_ng,gbleaf,gsleaf,co2ci)    ! output

    implicit none
!    Switches
    logical,parameter:: wrt_Aleaf=.false.
    
!    input variables
    logical wrt_photo_diag
    integer m,ng,yr,yr_spinup,phenoset
    real radabv(2)
    real Qabs(3,2)
    real Eabsorp_leaf_ng(2)
    real Cmolar,grn,psyc,slope
    real wind_l,Tair,TairK,Dair,co2ca,wleaf
    real raero,Ds0,swater_slw
    real Rconst,cpair,Patm,Temp_ref,H2OLv0,H2OLv,AirMa,H2OMw,Dheat
    real gsc0,alpha,stom_n
    real Vcmax_nl,Jmax_nl,Vcmax_par(3)
    real gam0,gam1,gam2
        
!   output variables      
    real Aleaf(2),Etransp_ng(2),Hleaf(2),Tleaf_ng(2),co2ci(2),gbleaf(2),gsleaf(2)      
        
!    internal variables
    real esat
    integer kr1,ileaf
    real Assimilation
    real gbw,gsc,gbc
    real rhocp,rrdn,TleafK_ng,TleafK_ng_tem
    real rbH,rbH_L,Y,co2cs,gbH,gbHf,gbHu,Grashof,gsw
    real Dleaf  ! Humidity deficit at leaf surface
    real Vcmax_nlt,Jmax_nlt
    
    
!    boundary layer conductance for heat - single sided, forced convection
!    (Monteith 1973, P106 & notes dated 23/12/94); E1 in Leuning et al. 1995

    gbHu=0.003*sqrt(wind_l/wleaf)    !m/s
    if(wind_l<=0.0)    gbHu=0.003 !*sqrt(-wind_l/wleaf)
  
    
    do ileaf=1,2              ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
        Tleaf_ng(ileaf)=Tair
        TleafK_ng=Tleaf_ng(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
        Dleaf=Dair*1000              !Pa        
!        first estimate for co2cs
        co2cs=co2ca               !mol/mol
        
!    ********************************************************************
        kr1=0                     !iteration counter for LE
        
        do   !iteration for leaf temperature, according to Leuning et al. 1995
            
            call VJmax_cal(Vcmax_nl,Jmax_nl,TleafK_ng,Rconst,Vcmax_par,    &   ! input
                &   Vcmax_nlt,Jmax_nlt)     ! output 
            
!            print*,"co2ci 3:",Qabs
                
            call photosynthesis(ileaf,phenoset,radabv,Qabs,Vcmax_nlt,Jmax_nlt, &
                &   CO2Cs,co2ci,Dleaf,TleafK_ng,Ds0,&
                &   swater_slw,gsc0,alpha,& ! input
                &   Rconst,Temp_ref,gam0,gam1,gam2,  & ! input
                &   Assimilation,gsc)  ! output
        
!           print*,"Aleaf(ileaf) 1:",ileaf,Aleaf(ileaf),Vcmax_nl,TleafK_ng
!            print*,"gsc:",gsc
            
!           Update leaf variables ------------------------------
!           Unit of assimilation: mol/m2/s
            Aleaf(ileaf) = Assimilation     
          
!           Calculate gbH ---------------------------------------------
!           gbH: unit of m/s; total boundary layer conductance to heat            
            Grashof=1.6e8*ABS(Tleaf_ng(ileaf)-Tair)*(wleaf**3)  !  Grashof number
!           Leaf boundary conductance of heat for free convection; E3 in Leuning et al. 1995 
            gbHf=0.5*Dheat*(Grashof**0.25)/wleaf
            gbH=gbHu+gbHf      
            gbw=1.075*gbH      !   Boundary layer conductance to water 
            
!           Calculate energy for transpiration ------------------------
!           Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            Y=1./(1.+ grN/gbH)
!           boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc=Cmolar*gbH/1.32  ! Calculate boundary layer conducatance of CO2 from boundary conductance of heat; mol/m2/s
            co2cs = co2ca-Aleaf(ileaf)/gbc
            co2Ci(ileaf)= co2cs-Aleaf(ileaf)/gsc
            if (co2Ci(ileaf).le.0) co2Ci(ileaf)= co2cs*0.10     ! to prevent a negative co2Ci(ileaf)

            gsw=gsc*1.56/Cmolar      ! stomatal conductance for H2O; unit of mol/m2/s

!           Etransp_ng: transpiration rate, according to Leuning 1995, Plant, Cell and Environment, 18: 1183-1200 
!           Dair: unit of Pa, according to Leuning 1995, Plant, Cell and Environment, 18: 1183-1200
            Etransp_ng(ileaf)= (slope*Y*Eabsorp_leaf_ng(ileaf)+rhocp*(Dair*1000)*gbH)/    &   !2* Weng 0215
            &     (slope*Y+psyc*gbH*(1/gbw+1/gsw)) 
!            print*,"Etransp_ng(ileaf):",
            
!          calculate sensible heat flux, according to Leuning 1995, Plant, Cell and Environment, 18: 1183-1200  
            Hleaf(ileaf)=Y*(Eabsorp_leaf_ng(ileaf)-Etransp_ng(ileaf))
!           gbH: leaf boundary-layer conductance to heat
            TleafK_ng_tem=273.2+Tair+Hleaf(ileaf)/(rhocp*gbH)
!           Dleaf: vapor deficit of leaf (Pa)
            Dleaf=psyc*Etransp_ng(ileaf)/(rhocp*gsw) 
            gbleaf(ileaf)=gbc*1.32*1.075
            gsleaf(ileaf)=gsc
    
!          compare current and previous leaf temperatures
            if(abs(TleafK_ng_tem-TleafK_ng).le.0.2) exit ! original is 0.05 C Weng 10/31/2008
!          update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
            TleafK_ng=TleafK_ng_tem
            Tleaf_ng(ileaf)=TleafK_ng_tem-273.2
            
!            print*,"Tleaf_ng:",ileaf,Etransp_ng(ileaf),Dair,gbH,psyc,gbH,gbw,gsw
            
            kr1=kr1+1
            if(kr1 > 100)then
                TleafK_ng=TairK
                exit
            endif
            if(TleafK_ng < 200.)then
                TleafK_ng=TairK
                exit 
            endif                     ! Weng 10/31/2008

        enddo
    enddo
!    print*,"Aleaf:",m,swater_slw,Aleaf(1),Aleaf(2)
    
    if(yr.gt.yr_spinup .and. wrt_aleaf)then
        write(3200,32001) m,ng,swater_slw,Aleaf(1),Aleaf(2)
!        print*,"Aleaf:",m,swater_slw,Aleaf(1),Aleaf(2)
    endif
32001   format(2(I8,","),2(f15.10,","),(f15.10))        
    
    return
end  ! leaf subroutine
!   *************************************************************************************


!****************************************************************************
subroutine photosynthesis(ileaf,phenoset,radabv,Qabs,Vcmax_nlt,Jmax_nlt, &
                    &   CO2Cs,co2ci,Dleaf,TleafK_ng,Ds0,&
                    &   swater_slw,gsc0,alpha,& ! input
                    &   Rconst,Temp_ref,gam0,gam1,gam2,  & ! input
                    &   Assimilation,gsc)  ! output
    implicit none
!    Switches
    logical,parameter :: photosyn_diagn=.false.  ! diag switch

    
!    input variable
    integer ileaf,phenoset
    real radabv(2)
    real Qabs(3,2)
    real CO2Cs ! leaf surface co2 concentration (mol/mol)
    real co2ci(2) ! leaf internal co2 concentration (mol/mol)
    real Vcmax_nlt,Jmax_nlt
    real Dleaf ! air vapour pressure deficit (KPa)  
    real TleafK_ng ! leaf temperature (oC)
    real Qapar ! quanta/m2/s
    real Ds0  ! Pa; An emperical constant reflecting the sensitivity of the stomata to Ds. Leuning et al. 1995
    real swater_slw ! fine root weighted water scaler
    real gsc0 ! residual stomatal conductance of co2 at the light compensation point
    real alpha
    real Kc_ref  ! Mechaelis constant for CO2 at constant temperature (25oC) 
    real Ko_ref  ! Mechaelis constant for O2 at constant temperature (25oC)
    real Ekc,Eko,o2ci,Rconst,Temp_ref,gam0,gam1,gam2

!   output variables
    real Assimilation   !
    real gsc     

!   internal variable
    real a1
    real Assim_rubi,Assim_rup2
    real Kc_temp,water_scal2assim
    real Rd,eJ,conKoT,gamma,gammas,tdiff
    real VJmax_temp,fJQres,VJtemp,VJtemp_Medlyn_peak
!   End of defining variables --------------------------------------------------------

    
!    print*,"co2ci:",ileaf,co2ci(ileaf)
    
!   Defining constants according to Leuning (1990) ------------
    Kc_ref=3.02E-04 ! Michaelis constant for CO2 at reference temperature (25oC); unit: mol/mol
    Ko_ref=2.56E-01 ! Michaelis constant for O2 at reference temperature (25oC); unit: mol/mol
    Ekc=5.94E+04    ! Activation energy for Kc; unit: J/mol
    Eko=3.60E+03    ! Activation energy for Ko; unit: J/mol
    o2ci=2.10E-01   ! O2 mol fraction; unit: mol/mol

    Qapar = (4.6e-6)*Qabs(1,ileaf)
    alpha=0.385	    ! unitless	A constant for calculating J (electron transfer rate)	

    
!   Calculate J, the asymptote for RuBP regeneration rate at given Q --------------
    eJ= fJQres(Jmax_nlt,alpha,Qapar)
    
!   Calculate Michaelis K ---------------------------------------------------------
!     Adjust Michaelis K for CO2 (Kc_temp) by temperature, according to Leuning 1990
    Kc_temp= Kc_ref*exp((Ekc/(Rconst* Temp_ref))*(1.-Temp_ref/TleafK_ng))

!     Adjust Michaelis K for CO2 (Kc_temp) by temperature
    conKoT=Ko_ref*exp((Eko/(Rconst* Temp_ref))*(1.-Temp_ref/TleafK_ng))

!     following de Pury 1994, eq 7, make light respiration a fixed proportion of
!   Rd: dark respiration    
    Rd= 0.0089*Vcmax_nlt                              !de Pury 1994, Eq7
!    print*,"Rd:",Rd,Vcmax_nlt 
    
    Tdiff=TleafK_ng-Temp_ref
!   gammas: co2 respiration point in the absence of mitochondrial respiration at 25oC??
!   calculate gammas according to Leuning et al. 1995.
!   gammas: CO2 compensation point without dark respiration; Farquhar et al, 1980; Luo and Reynolds, 1999       
    gammas= 34.6*(1.+0.0451*Tdiff+0.000347*Tdiff*Tdiff)/10e6       !gammas, unit of mol/mol
    
    gamma= (gammas+Kc_temp*(1.+O2ci/conKoT)*Rd/Vcmax_nlt)/(1.-Rd/Vcmax_nlt)
!    print*,"gamma:",gamma,gammas,Kc_temp,(1.+O2ci/conKoT),Rd,Vcmax_nlt,(1.-Rd/Vcmax_nlt)

!     calculate X using Lohammer model, and scale for soil moisture
!      range set according to Lohammer equation in Leuning 1995.
!      1/a1 = 1- ci/cs  
    a1=Amin1(10.0,Amax1(1.0,CO2Cs/(CO2Cs-co2ci(ileaf))))

!   Calculating stomatal conductance ! use water scaler
!      Ds0: an emperical coefficient reflecting the sensitivity of the stomata to Dleaf, unit of Pa
!      gamma: CO2 compensation point, unit of mol CO2/mol
!      a1: related to intercellular CO2 concentration at saturating irradiance
 
    water_scal2assim=Amin1(0.0+swater_slw,1.0)
        
!     calculate Rubisco activity limited assimilation (Assim_rubi)
!     according to equation C1 in Leuning et al. 1995; and Luo and Reynolds, 1999
    Assim_rubi=Vcmax_nlt*(co2ci(ileaf)-gammas)/(co2ci(ileaf)+Kc_temp*(1+o2ci/conKoT))
!    print*,"Assim_rubi:",Assim_rubi,Vcmax_nlt,co2ci(ileaf),gammas,Kc_temp,o2ci,conKoT 
    
!     calculate RuP2 regeneration limited Assim_rubi (Assim_rup2)
!     according to equation C2 in Leuning et al. 1995
!   Some difference with the equation (2) in Luo and Reynolds, 1999
    Assim_rup2=(eJ/4)*(co2ci(ileaf)-gammas)/(co2ci(ileaf)+2*gammas) 
    
    Assimilation= (amin1(Assim_rubi,Assim_rup2) - Rd)*water_scal2assim
    if(phenoset.eq.0) Assimilation=0.0    ! Assimilation*LAI

!    This works not bad: *Amin1(swater_slw+0.1,1.0) !

!    print*,"Assimilation:",Assimilation,Assim_rubi,Assim_rup2,Rd,eJ,Jmax_nlt,alpha,Qapar,Qabs(1,ileaf),ileaf
    
    if (isnan(Assimilation)) then
        print*,"Assimilation:",Assimilation,Assim_rubi,Assim_rup2,Rd
        print*,"Warning: Assimilation rate is NA!!!"
        stop
    endif
    
!    calculate new values for gsc (Lohammer model)
    gsc=gsc0+Assimilation*a1/((co2cs - gamma)*(1.0+Dleaf/Ds0)) !  ! revised by Weng
    gsc=Amax1(gsc,0.002)
!    print*,"Assimilation:",gsc,gsc0,X,Assimilation
    
    if(radabv(1).lt.10) then
        Assimilation=-Rd
        gsc=gsc0
    endif
    
    return
end  ! photosyn subroutine
!***********************************************************************

 !  Subroutine for calculating Vcmax_nlt and Jmax_nlt -------------------------------           
subroutine VJmax_cal(Vcmax_nl,Jmax_nl,TleafK_ng,Rconst,Vcmax_par,    &   ! input
                &   Vcmax_nlt,Jmax_nlt)     ! output 
    implicit none
!   input
    real Vcmax_nl,Jmax_nl,TleafK_ng,Rconst,Vcmax_par(3)
!   output
    real Vcmax_nlt,Jmax_nlt
!   internal
    integer methodtype
    real VJtemp_Medlyn_peak,VJmax_temp,VJtemp
    
!   Calculations --------------------------------------
    methodtype=1
    if(methodtype.eq.1)then ! Medlyn peak method
        Vcmax_nlt=VJtemp_Medlyn_peak(Vcmax_nl,TleafK_ng,Rconst,Vcmax_par,1)
        Jmax_nlt=VJtemp_Medlyn_peak(Jmax_nl,TleafK_ng,Rconst,Vcmax_par,2)
    else if(methodtype.eq.2) then   !  Leuning method
        Vcmax_nlt=VJmax_temp(Vcmax_nl,TleafK_ng,Rconst,1)
        Jmax_nlt=VJmax_temp(Jmax_nl,TleafK_ng,Rconst,2)
    else if(methodtype.eq.3)then    ! Reed optimal temperature method
        Vcmax_nlt=VJtemp(Vcmax_nl,TleafK_ng)
        Jmax_nlt=VJtemp(Jmax_nl,TleafK_ng)
    endif
    
    return
end
!   End of subtroutine VJmax_cal ===================================================




!****************************************************************************
real function esat(T)
    real T
!     returns saturation vapour pressure in Pa; Tetens equation
    esat=610.78*exp(17.27*T/(T+237.3))
    return
end

!****************************************************************************
      real function evapor(Td,Tw,Patm)
!* returns vapour pressure in Pa from wet & dry bulb temperatures
      gamma = (64.6+0.0625*Td)/1.e5
      evapor = esat(Tw)- gamma*(Td-Tw)*Patm
      return
      end

!****************************************************************************
      real function Vjmax(Tk,Temp_ref,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
      anum = Vjmax0*EXP((Eactiv/(Rconst*Temp_ref))*(1.-Temp_ref/Tk))
      aden = 1.+EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
      Vjmax = anum/aden
      return
      end
!****************************************************************************
      real function funE(extkbd,LAI)
      funE=(1.0-exp(-extkbd*LAI))/extkbd
      return
      end
!   ***************************************************************************
            
!     ****************************************************************************
!     Leuning et al. (1995) equation for temperature response
!     Medlyn et al. (2002)
!     used for Vcmax and Jmax
real function VJmax_temp(Vcmx1,TleafK_ng,Rconst,VJ)
!   input variable
    real Vcmx1,TleafK_ng,Rconst
    integer VJ
!   internal variable
    real Hv,Sv,Hd,T0
    real temp(3)
    
    if (VJ.eq.1) then       ! For Vcmax
        Hv=116300   ! unit of J/mol; initial: 116300
        Hd=202900   ! unit of J/mol; initial: 202900
    else if (VJ.eq.2) then  ! For Jmax
        Hv=79500   ! unit of J/mol
        Hd=201000   ! unit of J/mol; initial: 201000
    endif
    
    Sv=650  ! unit of J/mol; initial: 660
    T0=293.2    ! reference temperature; unit of K; initial: 293.2
    
    VJmax_temp=Vcmx1*exp((Hv/(Rconst*T0))*(1-T0/TleafK_ng))/(1+exp((Sv*TleafK_ng-Hd)/(Rconst*TleafK_ng)))
    
    return  
end
!   *****************************************************************  

real function VJtemp_Medlyn_peak(VJmax0,TleafK_ng,Rconst,Vcmax_par,VJ)
!   Method: peak function in Medlyn et al. (2002)
    implicit none
!   input
    real VJmax0,TleafK_ng,Rconst,Vcmax_par(3)
    integer VJ
    
!   internal
    real Ha,Hd,Topt
    real part1,part2
   
    Topt=Vcmax_par(1)
        
!    if (VJ.eq.1) then       ! For Vcmax
        Ha=Vcmax_par(2)   ! unit of J/mol; 
        Hd=Vcmax_par(3)   ! unit of J/mol;
!    else if (VJ.eq.2) then  ! For Jmax
!        Ha=60000   ! unit of J/mol
!        Hd=200000   ! unit of J/mol
!    endif

    
    part1=Hd*exp(Ha*(TleafK_ng-Topt)/(TleafK_ng*Rconst*Topt))
    part2=Hd-Ha*(1-exp(Hd*(TleafK_ng-Topt)/(TleafK_ng*Rconst*Topt)))
    
    VJtemp_Medlyn_peak=VJmax0*part1/part2

    return

end
!   ============================================================================

!     ****************************************************************************
!     Reed et al (1976, J appl Ecol 13:925) equation for temperature response
!     used for Vcmax and Jmax
real function VJtemp(VJmax0,TleafK_ng)
    implicit none
!   input
    real VJmax0,TleafK_ng
    
!   internal
    real TmaxVJ,TminVJ,ToptVJ,pwr
    
    
    TmaxVJ=50.0
    TminVJ=15.0
    ToptVJ=35.0
  
    if(TleafK_ng.lt.TminVJ) TleafK_ng=TminVJ   !constrain leaf temperatures between min and max
    if(TleafK_ng.gt.TmaxVJ) TleafK_ng=TmaxVJ
    pwr=(TmaxVJ-ToptVJ)/(ToptVj-TminVj)
    
    VJtemp=VJmax0*((TleafK_ng-TminVJ)/(ToptVJ-TminVJ))*     &
    &       ((TmaxVJ-TleafK_ng)/(TmaxVJ-ToptVJ))**pwr 
    return
end
!     End ===================================================================

!     ****************************************************************************
real function fJQres(Jmax_nlt,alpha,Qapar)
    real Jmax_nlt,alpha,Qapar,theta
    real AX,BX,CX
    
    theta=9.00E-01
    
    AX = theta                                 !a term in J fn
    BX = alpha*Qapar+Jmax_nlt                          !b term in J fn
    CX = alpha*Qapar*Jmax_nlt                          !c term in J fn
    if((BX*BX-4.*AX*CX)>=0.0)then
        fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
    else
        fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
    endif

    return
end
!     *************************************************************************


!     *************************************************************************
real function coszen_cal(doy,lat,hour)
    implicit none
!   inputs
    integer doy,hour
    real lat
!   internal variables
    real rad,pi,sinlat,coslat,sindec,cosdec,A,B
    
    
    !     sin(bet), bet = elevation angle  sun
    !     calculations according to Goudriaan & van Laar 1994 P30
    pi=3.1415926
    rad = pi/180.
    !     sine and cosine of latitude
    sinlat = sin(rad*lat)
    coslat = cos(rad*lat)
    !     sine of maximum declination
    sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
    cosdec=sqrt(1.-sindec*sindec)
    !     terms A & B in Eq 3.3
    A = sinlat*sindec
    B = coslat*cosdec
    coszen_cal = A+B*cos(pi*(hour-12.)/12.)
    return
end

!     *************************************************************************
subroutine fbeam_cal(doy,hour,lat,radsol,fbeam)
    implicit none
!   inputs
    integer doy,hour
    real lat,radsol
!   outputs
    real fbeam
!   internal variables
    real pi,pidiv,slatx,sindec,cosdec,a,b,sinbet,solext
    real tmprat,tmpR,tmpK,fdiff
    
    
    pi=3.14159256
    pidiv=pi/180.0
    slatx=lat*pidiv
    sindec=-sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
    cosdec=sqrt(1.-sindec*sindec)
    a=sin(slatx)*sindec
    b=cos(slatx)*cosdec
    sinbet=a+b*cos(2*pi*(hour-12.)/24.)
    solext=1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet

    tmprat=radsol/solext
    tmpR=0.847-1.61*sinbet+1.04*sinbet*sinbet
    tmpK=(1.47-tmpR)/1.66
   
    if(tmprat.le.0.22) fdiff=1.0
    if(tmprat.gt.0.22.and.tmprat.le.0.35) then
        fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
    endif
    if(tmprat.gt.0.35.and.tmprat.le.tmpK) then
        fdiff=1.47-1.66*tmprat
    endif
    if(tmprat.ge.tmpK) then
        fdiff=tmpR
    endif
    fbeam=1.0-fdiff
    if(fbeam.lt.0.0) fbeam=0.0
    return
end


!   Read parameters ========================================================================================
      
!   1. Read parameters and initial settings  ========================================================
      
!   1.1 Subroutine 1.1 Read parameters from file ============================================
subroutine Getpara(parafile,par_n,    &
            &   par_main,da_check,da_min,da_max)
    
    implicit none
!   input ---------
    integer par_n
    character(len=50) parafile
!   output ------------
    real par_main(par_n),da_min(par_n),da_max(par_n)
    integer da_check(par_n)
    character(len=50) commts
!    internal ----------
    integer m

    parafile=TRIM(parafile)
!   open and read input file for getting climate data
    
    open(1001,file=parafile,status='old')
    read(1001,10011)commts
10011  format(a132)

    do m=1,par_n   ! read parameters for C cycle
        read(1001,*) commts,commts,par_main(m),da_check(m),da_min(m),da_max(m)
    enddo ! end of reading the forcing file
    
    close(1001)
    return
end subroutine Getpara
!   1.1 End ==================================


!   ========================================================
!   Subroutine 1.2 Get parameters for soil_water subroutine
subroutine Getparameters_swater(parafile_swater,crt,wlama,potcof,wpot50,condb,n_pot)  ! parameter related to soil water redistribution
    
    implicit none
    real  crt,wlama,potcof,wpot50,condb,n_pot
    integer order
    character(len=50) parafile_swater,commts

    parafile_swater=TRIM(parafile_swater)
!   open and read input file for getting climate data
    
    open(1002,file=parafile_swater,status='old')
    read(1002,10021)commts
    read(1002,10021)commts
    read(1002,*)order,commts,crt
    read(1002,*)order,commts,wlama
    read(1002,*)order,commts,potcof
    read(1002,*)order,commts,wpot50
    read(1002,*)order,commts,condb
    read(1002,*)order,commts,n_pot

10021  format(a132)
    close(1002)
    return
    end  ! End of Getparameters_swater subroutine
!   1.2 End ==================================================
    
!   ========================================================
!   Subroutine 1.3 Get parameters for soil_temperature subroutine
    
subroutine Getparameters_stemp(parafile_stemp,shcap_snow,condu_snow,condu_b,&  ! parameter related to snow
    &   depth_ex,albedo_snow,resht, &   ! parameter related to snow
    &   fa,fsub,rho_snow,decay_m,shcap_soil,condu_soil,         &   ! parameter related to snow
    &   shcap_water,condu_water,shcap_air,condu_air,&
    &   shcap_ice,condu_ice,latent_heat_fusion,ice_density)  
    
    implicit none
    real  shcap_snow,condu_snow,condu_b,depth_ex,diff_s,albedo_snow,resht
    real  fa,fsub,rho_snow,decay_m
    real  shcap_soil,condu_soil,shcap_water,condu_water,shcap_air,condu_air
    real  shcap_ice,condu_ice,latent_heat_fusion,ice_density
    integer order
    character(len=50) parafile_stemp,commts

    parafile_stemp=TRIM(parafile_stemp)

    open(1003,file=parafile_stemp,status='old')
    read(1003,10031)commts
    read(1003,10031)commts
    read(1003,*)order,commts,shcap_snow              ! Order = 1
    read(1003,*)order,commts,condu_snow              ! Order = 2
    read(1003,*)order,commts,condu_b                 ! Order = 3
    read(1003,*)order,commts,depth_ex                ! Order = 4
    read(1003,*)order,commts,albedo_snow             ! Order = 7
    read(1003,*)order,commts,resht                   ! Order = 8
    read(1003,*)order,commts,fa                      ! Order = 12
    read(1003,*)order,commts,fsub                    ! Order = 13
    read(1003,*)order,commts,rho_snow                ! Order = 14
    read(1003,*)order,commts,decay_m                 ! Order = 15
    read(1003,*)order,commts,shcap_soil              ! Order = 16
    read(1003,*)order,commts,condu_soil              ! Order = 17
    read(1003,*)order,commts,shcap_water             ! Order = 18
    read(1003,*)order,commts,condu_water             ! Order = 19
    read(1003,*)order,commts,shcap_air               ! Order = 20
    read(1003,*)order,commts,condu_air               ! Order = 21
    read(1003,*)order,commts,shcap_ice               ! Order = 22
    read(1003,*)order,commts,condu_ice               ! Order = 23
    read(1003,*)order,commts,latent_heat_fusion      ! Order = 24
    read(1003,*)order,commts,ice_density             ! Order = 25

10031  format(a132)
    close(1003)
    return
end    ! End of Getparameters_stemp subroutine
!   1.3 End =====================================================
    
!   =============================================================
!   Subroutine 1.4 Get parameters for energy fluxes

subroutine Getparameters_energy(parafile_energy,tauL,rhoL,rhos,emleaf,emsoil,wleaf,gsc0)

    implicit none
    real,dimension(3):: tauL,rhoL,rhoS
    real emleaf,emsoil,wleaf,gsc0
    integer order
    character(len=51) parafile_energy,commts

    parafile_energy=TRIM(parafile_energy)
    
    open(1004,file=parafile_energy,status='old')
    read(1004,10041)commts
    read(1004,10041)commts
    read(1004,*)order,commts,tauL(1)
    read(1004,*)order,commts,rhoL(1)
    read(1004,*)order,commts,rhos(1)
    read(1004,*)order,commts,tauL(2)
    read(1004,*)order,commts,rhoL(2)
    read(1004,*)order,commts,rhos(2)
    read(1004,*)order,commts,tauL(3)
    read(1004,*)order,commts,rhoL(3)
    read(1004,*)order,commts,rhos(3)
    read(1004,*)order,commts,emleaf
    read(1004,*)order,commts,emsoil
    read(1004,*)order,commts,wleaf
    read(1004,*)order,commts,gsc0
    
10041  format(a132)
    close(1004)
    return
end     ! End of Getparameters_constants subroutine
!   1.4 End ====================================================



!   ========================================================
!   Subroutine 1.5 Get initial soil conditions 
    
subroutine Get_soilconditions(parafile_sconditions,  &
    &   layern,sthick,frlen,stemp,waterv,ice) 
    
    implicit none
    integer layern
    real sthick(10),frlen(10),stemp(10)
    real waterv(10),ice(10)
    integer order,i
    character(len=55) parafile_sconditions,commts

    parafile_sconditions=TRIM(parafile_sconditions)

    open(1005,file=parafile_sconditions,status='old')
    read(1005,10051)commts
    read(1005,10051)commts
    read(1005,10051)commts
    read(1005,10051)commts
    do i=1,layern
        read(1005,*)order,sthick(i),frlen(i),stemp(i),waterv(i),ice(i)
    enddo

10051  format(a132)
    close(1005)
    return
end     ! End of Getparameters_soilconditions subroutine 
!   1.5 End =================================================

!   ========================================================
!   Subroutine 1.6 Get initial pool related parameter values 
    
subroutine Getparameters_pools(parafile_pools,QC,CNini) 
    
    implicit none
    real  QC(8),CNini(8)
    integer order
    character(len=50) parafile_pools,commts

    parafile_pools=TRIM(parafile_pools)

    open(1006,file=parafile_pools,status='old')
    read(1006,10061)commts
    read(1006,10061)commts
    read(1006,10061)commts
    read(1006,10061)commts
    read(1006,*)order,commts,QC(1),CNini(1)      ! Order = 1
    read(1006,*)order,commts,QC(2),CNini(2)      ! Order = 2    
    read(1006,*)order,commts,QC(3),CNini(3)      ! Order = 3    
    read(1006,*)order,commts,QC(4),CNini(4)      ! Order = 4    
    read(1006,*)order,commts,QC(5),CNini(5)      ! Order = 5    
    read(1006,*)order,commts,QC(6),CNini(6)      ! Order = 6    
    read(1006,*)order,commts,QC(7),CNini(7)      ! Order = 7    
    read(1006,*)order,commts,QC(8),CNini(8)      ! Order = 8    

10061  format(a132)
    close(1006)
    return
end    !    End of Getparameters_pools subroutine 
!    1.6 End ==========================================
    
!   ========================================================
!   Subroutine 1.7 Get model initial values
    
subroutine Getparameters_siteconditions(parafile_siteconditions, &
    &   lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage)  
    
    implicit none
    real  lat,longi,Ttreat,CO2treat,Ndeposit,Nfert,Storage
    integer order
    character(len=55) parafile_siteconditions,commts

    parafile_siteconditions=TRIM(parafile_siteconditions)

    open(1007,file=parafile_siteconditions,status='old')
    read(1007,10071)commts
    read(1007,10071)commts
    read(1007,*)order,commts,lat                 ! Order = 1
    read(1007,*)order,commts,longi               ! Order = 2
    read(1007,*)order,commts,Ttreat              ! Order = 3
    read(1007,*)order,commts,CO2treat            ! Order = 4
    read(1007,*)order,commts,Ndeposit           ! Order = 5
    read(1007,*)order,commts,Nfert              ! Order = 6
    read(1007,*)order,commts,Storage             ! Order = 7
    
10071  format(a132)
    close(1007)
    return
end ! subroutine Getparameters
!   1.7 End ======================================================

    
!   =====================================================================================================
! Subroutine 2. Read climatic forcing from file   
subroutine Getclimate(climatefile,hour_length,clim_var_n, &  ! input
        &   year_seq,doy_seq,hour_seq,forcing_data,lines) ! output
    
    implicit none
    integer hour_length,clim_var_n
    integer,dimension(hour_length):: year_seq,doy_seq,hour_seq
    real forcing_data(clim_var_n,hour_length)
    character(len=150) climatefile,commts
    integer m,n,Reason,lines,year_length

    open(11,file=climatefile,status='old',ACTION='read',     &
    &     IOSTAT=Reason)
!    print*,"climatefile:",climatefile
    read(11,'(a160)') commts
!    print*,"test"
    
    m=0  ! to record the lines in a file
    year_length=0 ! to record years of a dataset
    
    do m=1,hour_length
!        m=m+1
        read(11,*,IOSTAT=Reason)year_seq(m),doy_seq(m),hour_seq(m),(forcing_data(n,m),n=1,clim_var_n)
!        if(Reason>0) then
!            print*,"Something is wrong!!!"
!        else if(Reason<0) then
!            print*,"End of the file reached"
!            exit
!        else
!        endif
    enddo

    lines=m-1
    year_length=(year_seq(lines)-year_seq(1))+1
!    print*,"Lines:",lines,year_length,forcing_data(1:1,1:10)
!    stop
    
    close(11)    ! close forcing file
    return
end
!   =====================================================================================================

! Subroutine 2. Read climatic forcing from file   
subroutine Getco2(co2file,co2_n,hour_length,co2data)
    
    implicit none
    integer hour_length,co2_n
    integer,dimension(hour_length):: year_seq,doy_seq,hour_seq
    real co2data(co2_n,hour_length)
    character(len=150) co2file,commts
    integer m,n,istat10,lines,year_length,test1,test2

    open(12,file=co2file,status='old',ACTION='read',     &
    &     IOSTAT=istat10)
    read(12,'(a160)') commts
    read(12,'(a160)') commts,test1,test2
!    print*,"commts:",commts,test1,test2
    m=0  ! to record the lines in a file
    year_length=0 ! to record years of a dataset
    do m=1,hour_length   ! read forcing files
!        m=m+1
        read(12,*,IOSTAT=istat10)year_seq(m),      &
        &       doy_seq(m),hour_seq(m),           &
        &       (co2data(n,m),n=1,co2_n),commts
        
!        print*,"CO2:",m,doy_seq(m),co2data(1,m),co2data(3,m),hour_length,co2file
    enddo ! end of reading the forcing file
    

    close(12)    ! close forcing file
    return
end
!   =====================================================================================================


    
!   =====================================================================================================
! Subroutine 3. read observation data from files

! read observed daily nee ==============================================
subroutine GetObs_d(obs_file,day_length,obs_d_n,obs_d_date_temp,obs_d_temp)
    implicit none
!    input
    character(len=150) obs_file
    integer day_length
    
!    output
    integer obs_d_date_temp(day_length),obs_d_n
    real obs_d_temp(day_length)
    
!    internal
    character(len=80) commts
    integer m,istat11
    
    
    open(2001,file=obs_file,status='old')  ! Open file
    
!    Read data
    read(2001,20011) commts
    obs_d_n=0
    do m=1,day_length
        read (2001,*,IOSTAT=istat11) obs_d_date_temp(m),obs_d_temp(m)  ! 1 days, 2 NEE
        if(istat11.lt.0)exit
        obs_d_n=obs_d_n+1
20011 format(a80)  
    enddo
    close(2001)

end subroutine GetObs_d 
! ===============================================================================


!   ==============================================================================
! Subroutine 1.1 Read estimated parameters      
Subroutine Getparaest(paraestfile,par_n,parapost,seq,npara,indexstring)
    implicit none
 
    character(len=50) paraestfile
    integer seq,m,n,istat6,par_n,npara
    real parapost(200,10000)
    character(len=par_n*15) indexstring,comments

    paraestfile=TRIM(paraestfile)
    open(15,file=paraestfile,status='old',ACTION='read',     &
    &     IOSTAT=istat6)

    read(15,*) comments
    read(15,'(A)') indexstring
    m=0
!   open and read input file for getting climate data
    do
        m=m+1
        read(15,*,IOSTAT=istat6)(parapost(n,m),n=1,npara)
        if(istat6<0)exit
    enddo
    seq=m-1

    close(15)
    return
end  !  subroutine Getparaest
!   ==============================================================================


!   ==============================================================================
subroutine getCov(gamma,covfile,npara)
    implicit none
    integer npara,i,k
    real gamma(npara,npara)    
    character(len=80) covfile
    
    open(14,file=covfile,status='old')

    do i=1,npara
        read (14,*)(gamma(i,k),k=1,npara)
    enddo  
    return
end ! subroutine getCov       
!   ==============================================================================
      
      
!   ==============================================================================
!    	cost function for observed data
    
subroutine Costfunction(isimu,Jscaler,day_length,day_length_ef,do_co2_da,do_swater_da, &  ! input
                &   obs_nee_d_n,obs_nee_d_date,obs_nee_d,simu_cflux_d,   & ! input
                &   obs_abc_d_n,obs_abc_d_date,obs_abc_d,   & ! input
                &   obs_swc2p5_d_n,obs_swc2p5_d_date,obs_swc2p5_d,simu_swc_d,   & ! input
                &   obs_swc12p5_d_n,obs_swc12p5_d_date,obs_swc12p5_d,   & ! input
                &   obs_swc22p5_d_n,obs_swc22p5_d_date,obs_swc22p5_d,   & ! input
                &   obs_swc37p5_d_n,obs_swc37p5_d_date,obs_swc37p5_d,   & ! input
                &   obs_swc52p5_d_n,obs_swc52p5_d_date,obs_swc52p5_d,   & ! input
                &   J_temp,J_last,upgraded,accR) ! output

    implicit none
!   Input
    logical do_co2_da,do_swater_da
    integer isimu,Jscaler,day_length,day_length_ef,obs_nee_d_n,obs_abc_d_n
    integer obs_swc2p5_d_n,obs_swc12p5_d_n,obs_swc22p5_d_n,obs_swc37p5_d_n,obs_swc52p5_d_n
    integer obs_nee_d_date(obs_nee_d_n),obs_abc_d_date(obs_abc_d_n),obs_swc2p5_d_date(obs_swc2p5_d_n)
    integer obs_swc12p5_d_date(obs_swc12p5_d_n),obs_swc22p5_d_date(obs_swc22p5_d_n)
    integer obs_swc37p5_d_date(obs_swc37p5_d_n),obs_swc52p5_d_date(obs_swc52p5_d_n)
    real obs_nee_d(obs_nee_d_n),obs_abc_d(obs_abc_d_n),simu_cflux_d(16,day_length_ef),simu_swc_d(5,day_length)
    real obs_swc2p5_d(obs_swc2p5_d_n),obs_swc12p5_d(obs_swc12p5_d_n),obs_swc22p5_d(obs_swc22p5_d_n)
    real obs_swc37p5_d(obs_swc37p5_d_n),obs_swc52p5_d(obs_swc52p5_d_n)
    
!   Output
    real J_temp,J_last,accR
    integer upgraded
    
!   Internal
    integer i,j,J_n
    real obs_nee_d_var,obs_abc_d_var,obs_swc_d_var(5)
    real r_num,J_new,delta_J
    real CalJ,CalVariance,J_nee_d,J_abc_d,J_swc_d(5)
    logical fixedvar,fixedprop
    
!   Initilization   
    fixedvar = .FALSE.
    fixedprop = .FALSE.
    if(fixedvar) then
        obs_nee_d_var= 0.5 ! g/m2/d
        obs_abc_d_var= 5  ! gC/m2
        obs_swc_d_var(1)=0.02  ! cm3/c3m3
        obs_swc_d_var(2)=0.02
        obs_swc_d_var(3)=0.02
        obs_swc_d_var(4)=0.02
        obs_swc_d_var(5)=0.02
    else if(fixedprop) then
        obs_nee_d_var= 0.0 ! g/m2/d
        obs_abc_d_var= 0.0  ! gC/m2
        obs_swc_d_var(1)=0.0  ! cm3/c3m3
        obs_swc_d_var(2)=0.0
        obs_swc_d_var(3)=0.0
        obs_swc_d_var(4)=0.0
        obs_swc_d_var(5)=0.0
    else
        obs_nee_d_var=CalVariance(obs_nee_d_n,obs_nee_d)
        obs_abc_d_var=CalVariance(obs_abc_d_n,obs_abc_d)
        obs_swc_d_var(1)=CalVariance(obs_swc2p5_d_n,obs_swc2p5_d)
        obs_swc_d_var(2)=CalVariance(obs_swc12p5_d_n,obs_swc12p5_d)
        obs_swc_d_var(3)=CalVariance(obs_swc22p5_d_n,obs_swc22p5_d)
        obs_swc_d_var(4)=CalVariance(obs_swc37p5_d_n,obs_swc37p5_d)
        obs_swc_d_var(5)=CalVariance(obs_swc52p5_d_n,obs_swc52p5_d)
    endif
    
    J_n=0
    J_nee_d=0.0
    J_abc_d=0.0
    J_swc_d=0.0

    if (do_co2_da) then
        if(obs_nee_d_n.gt.0) then
            J_nee_d=calJ(obs_nee_d_n,simu_cflux_d(1,obs_nee_d_date),obs_nee_d,obs_nee_d_var,fixedprop)
        endif
        if(obs_abc_d_n.gt.0) then
            J_abc_d=calJ(obs_abc_d_n,simu_cflux_d(4,obs_abc_d_date),obs_abc_d,obs_abc_d_var,fixedprop)
        endif
        J_n=J_n+1
    endif

    if (do_swater_da) then
        if(obs_swc2p5_d_n.gt.0) then
            J_swc_d(1)=calJ(obs_swc2p5_d_n,simu_swc_d(1,obs_swc2p5_d_date),obs_swc2p5_d,obs_swc_d_var(1),fixedprop)
        endif
        if(obs_swc12p5_d_n.gt.0) then 
            J_swc_d(2)=calJ(obs_swc12p5_d_n,simu_swc_d(2,obs_swc12p5_d_date),obs_swc12p5_d,obs_swc_d_var(2),fixedprop)
        endif
        if(obs_swc22p5_d_n.gt.0) then
            J_swc_d(3)=calJ(obs_swc22p5_d_n,simu_swc_d(3,obs_swc22p5_d_date),obs_swc22p5_d,obs_swc_d_var(3),fixedprop)
        endif
        if(obs_swc37p5_d_n.gt.0) then
            J_swc_d(4)=calJ(obs_swc37p5_d_n,simu_swc_d(4,obs_swc37p5_d_date),obs_swc37p5_d,obs_swc_d_var(4),fixedprop)
        endif
        if(obs_swc52p5_d_n.gt.0) then 
            J_swc_d(5)=calJ(obs_swc52p5_d_n,simu_swc_d(5,obs_swc52p5_d_date),obs_swc52p5_d,obs_swc_d_var(5),fixedprop)
        endif
        J_n=J_n+5
    endif

!    J_new=(J_nee_d+J_abc_d+sum(J_swc_d))/J_n  ! original method
    J_new=(J_nee_d+J_abc_d+sum(J_swc_d))/J_n*Jscaler
!    rescaler J_new value at upgraded 5/50/200 to obtain stable update rate, as controlled by Jscaler
!    if((upgraded.eq.50).or.(upgraded.eq.200) ) then
!        J_temp = Amax1(Amin1(J_temp*30/J_new,2.0),0.2)
!        J_new=(J_nee_d+J_abc_d+sum(J_swc_d))/J_n*30*J_temp
!        J_last = J_new
!    endif
    delta_J=J_new-J_last           !    delta_J=(J_new-J_last)/J_last
    
    CALL random_number(r_num)
    if(AMIN1(1.0,exp(-delta_J)).gt.r_num)then
        upgraded=upgraded+1
        J_last=J_new 
    endif
    
    accR=upgraded*1.0/isimu
    
    print*,isimu,upgraded,';delta_J:',int(delta_J*10)/10,";J new:",int(J_new*10)/10,";  acct rate(%):",int(accR*100)
            
    return
end   ! subroutine Costfunction
!   ==============================================================================


!   ==============================================================================
!   Square root of a matrix

subroutine racine_mat(M, Mrac,npara)

    integer npara,i
    real M(npara,npara),Mrac(npara,npara)
    real valpr(npara),vectpr(npara,npara)
    Mrac=0.
    call jacobi(M,npara,npara,valpr,vectpr,nrot)
    do i=1,npara
	if(valpr(i).ge.0.) then
            Mrac(i,i)=sqrt(valpr(i))
	else
            print*, 'WARNING!!! Square root of the matrix is undefined.'
            print*, ' A negative eigenvalue has been set to zero - results may be wrong'
            Mrac=M
            return
	endif
    enddo
    Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

end ! subroutine racine_mat      
!   ==============================================================================


!   ==============================================================================
!   Extraction of the eigenvalues and the eigenvectors of a matrix (Numerical Recipes)			
SUBROUTINE jacobi(a,n,np,d,v,nrot)
    INTEGER :: n,np,nrot
    REAL :: a(np,np),d(np),v(np,np)
    INTEGER, PARAMETER :: NMAX=500
    INTEGER :: i,ip,iq,j
    REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      
    do ip=1,n
        do iq=1,n
            v(ip,iq)=0.
        end do
        v(ip,ip)=1.
    end do

    do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
    end do

    nrot=0
    do i=1,50
        sm=0.
        do ip=1,n-1
            do iq=ip+1,n
                sm=sm+abs(a(ip,iq))
            end do
        end do
        if(sm.eq.0.)return
        if(i.lt.4)then
            tresh=0.2*sm/n**2
        else
            tresh=0.
        endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                    else
                        theta=0.5*h/a(ip,iq)
                        t=1./(abs(theta)+sqrt(1.+theta**2))
                        if(theta.lt.0.) then
                            t=-t
                        endif
                    endif
                    c=1./sqrt(1+t**2)
                    s=t*c
                    tau=s/(1.+c)
                    h=t*a(ip,iq)
                    z(ip)=z(ip)-h
                    z(iq)=z(iq)+h
                    d(ip)=d(ip)-h
                    d(iq)=d(iq)+h
                    a(ip,iq)=0.
                    do j=1,ip-1
                        g=a(j,ip)
                        h=a(j,iq)
                        a(j,ip)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                    end do
                    do j=ip+1,iq-1
                        g=a(ip,j)
                        h=a(j,iq)
                        a(ip,j)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                    end do
                    do j=iq+1,n
                        g=a(ip,j)
                        h=a(iq,j)
                        a(ip,j)=g-s*(h+g*tau)
                        a(iq,j)=h+s*(g-h*tau)
                    end do
                    do j=1,n
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
                    end do
                    nrot=nrot+1
                endif
            end do
        end do
        do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.
        end do
    end do
    print*, 'too many iterations in jacobi' 
    return
END  ! subroutine jacobi
!   ==============================================================================

!   ==============================================================================
!       generate new coefficents
subroutine gen_newcoef(coef_old,coefmax,coefmin,coef_new,search_length,npara)

    integer npara
    real coef_old(npara),coefmax(npara),coefmin(npara),coef_new(npara)
    real r,coefmid,random_harvest
    integer i
    real search_length
    do i=1,npara
999     continue
        CALL random_number(random_harvest)
        r=random_harvest-0.5
        coef_new(i)=coef_old(i)+r*(coefmax(i)-coefmin(i))*search_length
        if(coef_new(i).gt.coefmax(i).or.coef_new(i).lt.coefmin(i))goto 999
    enddo
    return
end  ! subroutine gen_newcoef
!   ============================================================================== 

!   ==============================================================================
!   Generation of a random vector from a multivariate normaldistribution 
!   with mean zero and covariance matrix gamma.									  !!
!   Beware!!! In order to improve the speed of the algorithms, 
!   the subroutine use the Square root matrix of gamma.

subroutine gengaussvect(gamma_racine,xold,xnew,npara)

    integer npara
    real gamma_racine(npara,npara)
    real x(npara),xold(npara),xnew(npara)

    do i=1,npara
        x(i)=rangauss(25)
    enddo

    x = matmul(gamma_racine, x)
    xnew = xold+x
end ! subroutine gengaussvect
!   ==============================================================================

!   ==============================================================================
!   Generation of a random number from a standard normal distribution. (Numerical Recipes)           !!

function rangauss(idum)


    integer idum
    real v1, v2, r, fac, gset
    real r_num

    data iset/0/
    if(iset==0) then
    1	CALL random_number(r_num)
            v1=2.*r_num-1
            CALL random_number(r_num)
            v2=2.*r_num-1
            r=(v1)**2+(v2)**2
            if(r>=1) go to 1
            fac=sqrt(-2.*log(r)/r)
            gset=v1*fac
            rangauss=v2*fac
            iset=1
    !
    else
        rangauss=gset
            iset=0
    end if

    return
end 
!   ==============================================================================


!   ==============================================================================
!! Compute the centered matrix, ie. the matrix minus the column means									  !!

subroutine centre(mat,mat_out,npara,ncov)
    integer npara,ncov
    real mat(ncov,npara),mat_out(ncov,npara)
    real mean
    do i=1,npara
        mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
    enddo
end !  subroutine centre
!   ==============================================================================


!   ==============================================================================
!! mean of a vector									  !!

Function mean(tab,ncov)
    integer ncov
    real tab(ncov)
    real mean,mean_tt
    mean_tt=0.
    do i=1,ncov	
        mean_tt=mean_tt+tab(i)/real(ncov)
    enddo
    mean=mean_tt
End Function

!  ===============================================================================
! Calculate the J value of each measurement

Real Function  CalJ(Size,Simulation,Observation,Variance,fixedprop)
    Implicit None
    Integer Size,i
    Real Simulation(size),Observation(size)
    Real Variance, Jtemp, Varprop
    Logical fixedprop
    
!    variation as a proportion of observed value (Varprop)
    Varprop = 0.5
    Jtemp = 0.0
    if(fixedprop) then
        do i=1,Size
            Jtemp=Jtemp+(Simulation(i)-Observation(i))**2/(Observation(i)*Varprop*1)
        enddo
        CalJ = Jtemp/Size
    else
!        CalJ=(NORM2(Simulation-Observation))**2/(Variance)  ! original
        CalJ=(NORM2(Simulation-Observation))**2/(Variance*Size)  ! standardized by sample size
    endif
    
END Function

!  ===============================================================================
! Calculate the standard deviation of an array

Real Function  CalVariance(Size,Dataset)
    Implicit None
    Integer Size, i, size2,size3
    Real Dataset(Size)
    Real Var_acc, Mean, Mean_acc

    Mean_acc = 0.0
    size2=0
    DO i = 1, Size
        if (Dataset(i).gt.-999) then
            Mean_acc = Mean_acc+Dataset(i)
            size2=size2+1
        endif
    END DO
    Mean = Mean_acc/Size2

    Var_acc = 0.0
    size3=0
    DO i = 1, Size
        if (Dataset(i).gt.-999) then
            Var_acc = Var_acc+(Dataset(i) - Mean)**2
            size3=size3+1
        endif
    END DO
    CalVariance = Var_acc /(size3-1)
    
END Function

