!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module usersettings
implicit none


    ! Define parameters for the model grid

    integer, parameter :: nx = 600! Number of x points (orig. 83) !need nx=125 for dom width 50km
    integer, parameter :: nz = 170 ! Number of z points orig. 42
    integer, parameter :: nt = 10000 !number of timesteps we want to run (for dt =2.0, need nt=600 to get to 1200s, for dt=0.5, need nt=2400)
    real, parameter :: dt = 0.5 !length of timestep in seconds (orig 2)
    real, parameter :: dx = 200.0 ! Horizontal grid spacing (m) 
    real, parameter :: dz = 100.0 ! Vertical grid spacing (m)
    real, parameter :: dtout = 60 ! Time spacing of output (used in writing out) (orig. 2)
	integer, parameter :: ntout_basest = 1 !number of times written out for base state vars (only one bc they dont change)
    integer, parameter :: ntout_perts = 10001! # of times written out for pert vars (will be calculated later)(eg. need to have ntout_perts = 601 if want to run 1200s with 2.0s timstep)

    real, parameter :: cs = 50.0 !speed of sound (m/s)
	! define variables/values for warm bubble
    real, parameter :: xc = 30000.0 !(meters, x coord of bubble's center point)(16200m is the domain's center scalar point)
    real, parameter :: zc = 3000.0 !(meters, z position of bubble center pt) orig 3000
    real, parameter :: rx = 4000.0 !(meters, bubble's half-width/radius) orig. 4000
    real, parameter :: rz = 4000.0 !(meters, bubble half-depth) orig. 4000
    real, parameter :: delT = 3.0 !(Kelvin, bubble's maximum temperature perturbation) !orig. 3.0 for warmbubble

        ! Define format specifier (specify how data is formatted in output)
    character(len=80) :: FMT

    !specify unitless cmixh and cmixv, used later to compute kmixh and kmixv 
    !originally used 0.01
    real, parameter :: cmixh = 0.01
    real, parameter :: cmixv = 0.01

    !specify the magnitude we want our random noise to be
    real, parameter :: amplitude = 0.0

    !add variables to specify for model top rayleigh damping
    real, parameter :: raydmpc = 0.01 !desired max damping rate (ar in the eqn)
    real, parameter :: zbottom = 11000 !height of the bottom damping layer (in m)

    !specify time filtering coefficient
    real, parameter :: epsilon = 0.1

    !specify parameters for surface fluxes
    real, parameter :: fluxmag = 0.002 ! average surface flux in K/s
    real, parameter :: snoise = 0.5 !random surface noise, so that turbulent motions will arise
    !snoise is used as a unitless fraction of fluxmag (so =0.5 means that noisy surface fluxes would range 50% above and below fluxmag)

    !SWITCH FOR SURFACE FLUXES (1 = on, 0 = off)
    real, parameter :: surface_fluxes = 0.0

    !TERRAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !specify bounds of terrain (USING GRIDPOINTS)
    !ex. if maxz_terr = 12.0, then the terrain ends at the 12th vertical gridpt (the height of that depends on grid spacing)
    !for simple 'plateau mtn'
    real, parameter :: minx_terr = 90.0 !x values that will be terrain (masked as =11)
    real, parameter :: maxx_terr = 100.0
    
    real, parameter :: minz_terr = 1.0 !this shouldn't change.... unless want a floating mtn lol
    real, parameter :: maxz_terr = 12.0
    
    !VARIABLES FOR STEP MTN (using formula from Rotunno and Bryan 2018)
    !need to be in meters, so multiply by the grid spacing
    !SWITCH FOR IF WANT TO USE MULTIPLE MTNS Or single mtn
    
    real, parameter :: multi_mtn = 0.0 ! 1 = multiple mtns, 0= single mtn
    
    !switch above determines if the '__2, etc.. parameters will be used in terrain_setup
    real :: a2 = 0.0
    real :: x_cent2 = 0.0
    real :: hmax2 = 0.0

    real :: a3 = 0.0
    real :: x_cent3 = 0.0
    real :: hmax3 = 0.0

    real :: a4 = 0.0
    real :: x_cent4 = 0.0
    real :: hmax4 = 0.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !SINGLE MTN
    real, parameter :: a = 50.0 * dx !half width of the mtn in gridpoints (so total mtn will be a*2 gridpoints in width) !orig 40
    real, parameter :: x_cent = 300.0 * dx!center gridpoint of the mtn (orig. 300)
    real, parameter :: hmax = 20.0 * dz !maximum height of the mtn in number of gridpoints (for reference, mt mitchell is ~2000m, so 20)

    ! real, parameter :: a = 0.0 * dx !half width of the mtn in gridpoints (so total mtn will be a*2 gridpoints in width)
    ! real, parameter :: x_cent = 0.0 * dx!center gridpoint of the mtn
    ! real, parameter :: hmax = 0.0 * dz !maximum height of the m

end module usersettings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module declarations

	use usersettings
	implicit none

    !first declare d2t, this is the variable we will use in the leapfrog scheme for pert vars(since leaprog scheme is 2*dt)(updates automatically when adjust dt)

    real :: d2t = 1*dt !set as 1dt until the end of the first timestep (at the last time in time loop, set this var to 2dt)(the first timestep is only forward, so need 1dt)
    
	! Physical constants

    real, parameter :: g = 9.8 ! Gravitational acceleration (m/s^2)
    real, parameter :: cp = 1004.0 ! Specific heat at constant pressure (J/kgK)
    real, parameter :: rd = 287.0 ! Dry gas constant (J/kgK)
    real, parameter :: cv = cp - rd ! Specific heat at constant volume (J/kgK)
    real, parameter :: trigpi = 4.0 * atan(1.0) ! Trigonometric pi
    real, parameter :: lv = 2.5E6 !latent heat of vaporization (J/kg)

    real, parameter :: p0 = 100000 !p0 in Pa (for use in calculations) 

    ! Declare the variables used to define grid point locations
    ! Initialize each array as zero for safety/debugging
    ! The grid setup/stagger is such that the w and scalar points have the same x-location,
    ! and the u and scalar points have the same z-location

    real, dimension(nx) :: xu = 0.0 ! x-location of u-points (m)
    real, dimension(nx) :: xs = 0.0 ! x-location of w-points & scalar points (m)
    real, dimension(nz) :: zw = 0.0 ! Height/z-position of w-points (m)
    real, dimension(nz) :: zs = 0.0 ! Height of u-points & scalar points (m)

    ! Define r for the warm bubble
    real, dimension(nx, nz) :: r = 0.0


	! Variable declarations (these are base state variables and are only a function of height)

    real, dimension(nz) :: ub = 0.0 ! Mean u-wind in m/s (ubar)
    real, dimension(nz) :: thb = 0.0 ! Mean potential temperature in K (thetabar)
    real, dimension(nz) :: pib = 0.0 ! Mean non-dimensional pressure (unitless)
    real, dimension(nz) :: qb = 0.0 ! Mean mixing ratio (kg/kg)
    real, dimension(nz) :: thvb = 0.0 ! Mean virtual potential temperature (K)
    real, dimension(nz) :: rhos = 0.0 ! Density on u-points & scalar points (kg/m^3) (u&s have the same vertical coordinate)
    real, dimension(nz) :: rhow = 0.0 ! Density on w-points (kg/m^3)
    real, dimension(nz) :: p = 0.0 ! Physical base pressure in Pa (can convert the non-dim pressure into this)
    !create a variable that contains values of thb but in an array filled for all times
    !use this to calculate total theta (base state + perts)
    real, dimension(nz) :: thb_alltime = 0.0

    integer :: i
    integer :: k
    integer :: n
    
    !declare total virtual potential temp (full thetav)(will use to calculate thetav')
    real, dimension(nx,nz) :: thvtot = 0.0 !full thetav (thvb + thv)

	! Declare perturbation variables (need 3 distinct 'time levels' for each pert variable, to use leapfrog differencing scheme)
    real, dimension(nx, nz) :: thm = 0.0 ! Pot temp. pert at t-dt (m/s) (m for minus)
    real, dimension(nx, nz) :: th = 0.0 ! Pot temp pert at t (m/s)
    real, dimension(nx, nz) :: thp = 0.0 ! Pot temp pert at t+dt (m/s) (p for plus)

    real, dimension(nx, nz) :: um = 0.0 ! U-wind pert at t-dt
    real, dimension(nx, nz) :: u = 0.0 ! U-wind pert at t
    real, dimension(nx, nz) :: up = 0.0 ! U-wind pert at t+dt

    real, dimension(nx, nz) :: wm = 0.0 ! Vertical wind pert at t-dt
    real, dimension(nx, nz) :: w = 0.0
    real, dimension(nx, nz) :: wp = 0.0

    real, dimension(nx, nz) :: pim = 0.0 ! Non-dimensional pressure pert at t-dt
    real, dimension(nx, nz) :: pi = 0.0
    real, dimension(nx, nz) :: pip = 0.0

    real, dimension(nx,nz) :: ppert = 0.0 !dimensional (physical) pert pressure
    
    real, dimension(nx,nz) :: pitot = 0.0 ! total non dim pressure

    real, dimension(nx,nz) :: ptot = 0.0 !total dim (physical) pressure

    real, dimension(nx, nz) :: qvm = 0.0 ! Water vapor mixing ratio perturbation at t-dt
    real, dimension(nx, nz) :: qv = 0.0
    real, dimension(nx, nz) :: qvp = 0.0

    real, dimension(nx, nz) :: qcm = 0.0 !condensate mixing ratio
    real, dimension(nx, nz) :: qc = 0.0
    real, dimension(nx, nz) :: qcp = 0.0

    real, dimension(nx, nz) :: psim = 0.0 !tracer that is passively transported by the wind
    real, dimension(nx, nz) :: psi = 0.0
    real, dimension(nx, nz) :: psip = 0.0 

    real, dimension(nx,nz) :: thvm = 0.0 !virtual potential temp 
    real, dimension(nx,nz) :: thv = 0.0 !virtual potential temp at present time
    real, dimension(nx,nz) :: thvp = 0.0

    !total theta (basestate plus pert)
    real, dimension(nx,nz) :: thtot = 0.0

	!define vars for the u and v inetrpolated to scalar pts
	real, dimension(nx,nz) :: uinterp = 0.0
	real, dimension(nx,nz) :: winterp = 0.0

    real, dimension(nx,nz) :: up_interp = 0.0 !future u values/precidted values interpolated
    real, dimension(nx,nz) :: wp_interp = 0.0 !predicted w values interpolated

    real :: pisfc ! Surface pressure (non-dimensional pressure at the physical surface) (will use to calculate pib later)

    real :: thvbavg ! Average thetav for each layer (will use later in pib calculation)

    ! Declare variables that will be used when reading in our input sounding

    real :: weight ! Use this for interpolating sounding data onto model grid points

    real :: inp ! Stores pressure
    integer, parameter :: num_z = 1000 ! Number of levels to read in from sounding
    real, dimension(num_z) :: inz, inth, inqv, inu, inzagl
    integer :: kk, nlevs ! kk is the counter for levels read in, nlevs stores the number of levels read

    !declare mixing coefficients for spatial filtering/diffusion

    real :: kmixh = 0.00
    real :: kmixv = 0.00

    !declare varibale to store random number from the random number generator
    real :: rand = 0.0

    !define variable for damping coefficient (tau in the equation)
    real :: coef = 0.0

    !define variables that will use for saturation adjustment (produce condensate when supersaturation occurs)
    
    real :: pres_tot = 0.0 !total physical pressure, calculated with future values (for instantaneous saturation adjustment)
    real :: temp_tot = 0.0 !total temp, calculated w future values
    real :: qv_tot = 0.0 !total vapor pressure, calculated w future values

    real :: qvsat = 0.0 !saturation vapor pressure (will be computed using Teten's formula)
    real :: phi = 0.0 !adjustment factor (accounts for changes in qvsat that occur as latent heating/cooling occurs (& thus temp changes))
    real :: dqv = 0.0 ! change in total water vapor due to evap or condensation (computed later using qvsat, phi, etc.)

    !variables for testing qc and qv
    real :: qv_one = 0.0
    real :: qv_two = 0.0
    real :: qv_three = 0.0

    real :: qc_one = 0.0
    real :: qc_two = 0.0

    real :: th_one = 0.0
    real :: th_two = 0.0
    real :: th_three = 0.0

    !declare variables needed for fallout of precip
    real, dimension(nx,nz) :: vt = 0.0 !terminal velocity of condensate fallout

    !declare vars for adding terrain
    real, dimension(nx,nz) :: terrain =0.0

    real, dimension(nx,nz) :: terrain_mask = 0.0

    real, dimension(nx,nz) :: mtnsurf = 0.0
    !multiple mtns
    real, dimension(nx,nz) :: mtnsurf2 = 0.0 

    real, dimension(nx,nz) :: mtn = 0.0
    !multiple mtns
    real, dimension(nx,nz) :: mtn2 = 0.0

    real :: elevation1 = 0.0
    real :: elevation2 = 0.0
    real :: elevation3 = 0.0
    real :: elevation4 = 0.0
    !DECLARE FUNCTION FOR CALCULATING TERRAIN
    contains

    !real function calculate_mountain_surface(xs, x_cent, a, hmax)
    real function calculate_mountain_surface(scalarpt, mtncent, hlfwidth, maxh)
        real, intent(in) :: scalarpt, mtncent, hlfwidth, maxh
    
        mtn(i,k) = ABS(scalarpt - mtncent) / hlfwidth
        
        if (mtn(i,k) .le. 1.0) then
            calculate_mountain_surface = maxh * cos(trigpi/2 * ABS(scalarpt - mtncent)/hlfwidth)**2
        else
            calculate_mountain_surface = 0.0
        endif
    end function calculate_mountain_surface
    
    
    
end module declarations
program packdriver

use declarations

implicit none
!setting up grid, ICs, opening output files, writing headers and base state vars
!only need to do this once at beginning of run
call setup_outfiles
call physical_grid_locs
call terrain_setup
call read_sounding
call compute_basestate

call terrainBCs !switch/ if check for terrain inside subroutine
call writeout_basestate
!call warmbubble
call initial_rand_noise
call leapfrog_setup
call zerogradLB_pert_var_BCs
!call periodicLB_pert_var_BCs
write(fmt, '(a,i4,a)') '(', nx, '(e10.4,1x))' ! Write out all of the x values on one line first
call writeout_perts

do n=1,nt 
    print *, '-----------------------------------------------------------'
    print *, 'STARTING TIMESTEP', n
    !for the first step, we want to only do a step forward of 1*dt (this is because our first timestep is technically just a forward step (due to the initial conditions set above))
    !it is technically just a forward step because since at the first step we are setting past=present, we are really just going from present to future instead of past to future)
    call flux_predictive_eqns

    !zero gradient in z when doing vertical diffutuon (do it before vert diff)
    !same thing for horiz
    call terrainBCs

    call spatial_filtering !includes the zero gradient along terrain faces for diffusion... if check

    call terrainBCs

    call saturation_adjust
    call modeltop_damping
     
    !Boundary conditions.. one subroutine has periodic LBCs and one has zero gradient (select which to run)
    call zerogradLB_pert_var_BCs
    !call periodicLB_pert_var_BCs
    call timefilter
    call handoff_values
    call writeout_perts
    !we want to make sure d2t is set to 2*dt, which we will need for the leapfrog scheme for all timesteps after the first timestep
    !resets the timestep to 2*dt after the first timestep and after every subsequent one
    d2t = 2*dt !now the timestep changes to 2dt for all timesteps after the very first one
    print *, 'FINISHED timestep', n
    print *, '-----------------------------------------------------------'
enddo

call close_files

stop
endprogram packdriver!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_outfiles
    use declarations
    implicit none

    ! open text files to store varaible data so that they can be plotted, and set up headers for outfiles
    ! Write out a different text file for each variable

    open(10, file='ub_out.txt')    ! ub 
    write(10, *) '# NX DX NZ DZ NTOUT DTOUT'   ! ub 
    write(10, *) '# ', nx, dx, nz, dz, ntout_basest, dtout     ! ub

    open(11, file='thb_out.txt')     
    write(11, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(11, *) '# ', nx, dx, nz, dz, ntout_basest, dtout     

    open(12, file='pib_out.txt')     
    write(12, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(12, *) '# ', nx, dx, nz, dz, ntout_basest, dtout

    open(13, file='qb_out.txt')     
    write(13, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(13, *) '# ', nx, dx, nz, dz, ntout_basest, dtout

    open(14, file='th_out.txt')     
    write(14, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(14, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(15, file='u_out.txt')     
    write(15, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(15, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(16, file='w_out.txt')     
    write(16, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(16, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(17, file='ppert_out.txt')     
    write(17, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(17, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(18, file='qv_out.txt')     
    write(18, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(18, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(19, file='qc_out.txt')     
    write(19, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(19, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

	open(20, file='p_out.txt') !dimensional pressure
	write(20, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(20, *) '# ', nx, dx, nz, dz, ntout_basest, dtout

    open(21, file='thv_out.txt') !theta v prime
	write(21, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(21, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(22, file='pi_out.txt') 
	write(22, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(22, *) '# ', nx, dx, nz, dz, ntout_basest, dtout

    open(23, file='psi_tracer_out.txt')     
    write(23, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(23, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(24, file='thtot_out.txt')     
    write(24, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(24, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    open(25, file='qv_tot_out.txt')     
    write(25, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(25, *) '# ', nx, dx, nz, dz, ntout_perts, dtout


    open(26, file='terrain_mask.txt')     
    write(26, *) '# NX DX NZ DZ NTOUT DTOUT'    
    write(26, *) '# ', nx, dx, nz, dz, ntout_perts, dtout

    

return
end subroutine setup_outfiles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine terrain_setup
    use declarations
    implicit none

    !create a simple terrain mask (plateau mtn/ tophat mtn... step mtn below)
    ! do i=1, nx
    !     do k=1, nz
    !         if ((i .ge. minx_terr) .and. (i .le. maxx_terr) .and. (k .ge. minz_terr) .and. (k .le. maxz_terr)) then
    !             terrain(i,k) = 1.0
    !         else
    !             terrain(i,k) = 0.0
    !         endif
    !     enddo
    ! enddo
    
    ! includes code for multiple mtns, with switch for single mtn if desired only one
    ! switch is in usersettings
    if (multi_mtn == 1.0) then
        !adjust these based on where want 2nd mtn to be
        a2 = 36.0 * dx 
        x_cent2 = 320.0 * dx !orig 220
        hmax2 = 16.0 * dz !imitate peak to west of richland basalm (on other side of that valley) (~1800m)

        a3 = 40 * dx
        x_cent3 = 390 * dx
        hmax3 = 16.5 * dz !use 18 for imitating richland basalm/ easternmost mtn of great basalm mtns (~1900m)

        a4 = 40.0 * dx
        x_cent4 = 290.0 * dx
        hmax4 = 16.7 * dz
    else
        a2 = 0.0
        x_cent2 = 0.0
        hmax2 = 0.0

        a3 = 0.0
        x_cent3 = 0.0
        hmax3 = 0.0

        a4 = 0.0
        x_cent4 = 0.0
        hmax4 = 0.0
    endif
    
    ! the function (calculate_mountain_surface)is defined within module (cannot be defined within subroutine)
    do i=1,nx
        do k=1, nz
            ! Calculate elevation for each mountain
            elevation1 = calculate_mountain_surface(xs(i), x_cent, a, hmax) !orig mtn(domain center, change height in usersettings)
            elevation2 = calculate_mountain_surface(xs(i), x_cent2, a2, hmax2) !second from left mtn
            elevation3 = calculate_mountain_surface(xs(i), x_cent3, a3, hmax3) !rightmost mtn
            elevation4 = calculate_mountain_surface(xs(i), x_cent4, a4, hmax4) !leftmost mtn
    
            ! Combine elevations of all terr features
            terrain(i,k) = max(elevation1, elevation2, elevation3, elevation4)
    
            ! Set the terrain mask based on combined elevation of all mtns
            if (zs(k) .le. terrain(i,k)) then
                terrain_mask(i,k) = 1
            else
                terrain_mask(i,k) = 0
            endif
        enddo
    enddo
    
    ! do i=1,nx
    !     do k=1, nz
    !         mtn(i,k) = ABS(xs(i) - x_cent) / a
    !         !mtn2(i,k) = ABS(xs(i) - x_cent2) / a2

    !         if (mtn(i,k) .le. 1.0) then
    !             mtnsurf(i,k) = hmax * cos(trigpi/2 * ABS(xs(i) - x_cent)/a)**2
    !             !mtnsurf2(i,k) = hmax2 * cos(trigpi/2 * ABS(xs(i) - x_cent2)/a2)**2
    !         else
    !             mtnsurf(i,k) = 0.0
    !         endif
    !     enddo
    ! enddo

    ! do i=1, nx
    !     do k=1, nz
    !         if (zs(k) .le. mtnsurf(i,k)) then !if zscalar is below (less than) the calculated mtn surface, then the point is terrain (underground)
    !             terrain(i,k) = 1
    !         else
    !             terrain(i,k) = 0
    !         endif
    !     enddo
    ! enddo
    ! !after predict vars... if terrain =1 then u, w, s =0

return
end subroutine terrain_setup!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine terrainBCs
    use declarations
    implicit none


    !set the initial fields to zero before the time loop
    !set all predicted variables to zero on all terrain faces (no values 'in the ground')
    !do the 'plus'(future) values and current values to be safe

    !set all 4 'faces' to zero (with respect to a given i,j point which is classified as terrain)
    !points above and to the right of the ground should be zero
    do i=2, nx-1
        do k=2, nz-1
            if (terrain_mask(i,k) .ge. 1.0) then
                u(i+1,k) = 0 
                up(i+1,k) = 0
                u(i-1,k) = 0
                up(i-1,k) = 0
                u(i,k+1) = 0
                up(i,k+1) = 0
                !u(i,k-1) = 0
                !up(i, k-1) = 0 !do an if check? if k is outside phys domain then dont do
            
                w(i+1,k) = 0 
                wp(i+1,k) = 0
                w(i-1,k) = 0
                wp(i-1,k) = 0
                w(i,k+1) = 0
                wp(i,k+1) = 0
                !w(i,k-1) = 0
                !wp(i, k-1) = 0
            endif
        enddo
    enddo
 
    !since w and u are zero on the terrain faces, we just need to keep setting the scalars inside the mtn to zero 
    !w and u being zero means nothing gets fluxed in or out of the mtn

    !set underground values of scalars = 0
    do i =2, nx-1
        do k=2, nz-1
            if (terrain_mask(i,k) .ge. 1.0) then
                thp(i,k) = 0
                th(i,k) = 0
                qv(i,k) = 0
                qvp(i,k) = 0
                pi(i,k) = 0
                pip(i,k) = 0
                psi(i,k) = 0
                psip(i,k) = 0
                qc(i,k) = 0
                qcp(i,k) = 0
            endif
        enddo
    enddo

return
end subroutine terrainBCs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine terrain_filt_horiz(varm)
    use declarations
    implicit none

    !modify spatial filtering to prevent things from diffusing in/out of terrain face/ground
    !need this so that points inside the terrain are not used for diffusion (or in the centered diff)
    !prevents loss of signal into mtn
    !checks when crossing into/out of the mtn (if we are at a physical point that is next to a terrain pt)
    !take that artificial point in the mtn, and make it equal to the last physical point
    
    !then do spatial filtering (as it is) -> (zero gradient must be done before the horizontal diffusion, so that things arent diffused into the mtn)
    
    real, dimension(nx, nz) :: varm
    !for a point i,j in the model atmopshere... 

    if (terrain(i,k) .eq. 0.0 .and. terrain(i+1,k) .eq. 1.0) then !if the point to the right of a point is inside the terrain mask... apply zero gradient condition
        varm(i+1,k) = varm(i,k) !this is when there is terrain to the left of point i,j.. set the point inside the mtn to the same value as the last physical point
        !we do this so that the point with 'zero' value inside the mtn does not factor into diffusion and decrease values of physical pts
    endif

    if (terrain(i,k) .eq. 0.0 .and. terrain(i-1,k) .eq. 1.0) then !if the point to the left of a point is inside the terrain mask... apply zero gradient condition 
        varm(i-1,k) = varm(i,k)
    endif

            !address possibility that mtn top could be 1 grid cell wide
    if (terrain(i,k) .eq. 1.0 .and. terrain(i-1,k) .eq. 0.0 .and. terrain(i+1,k) .eq. 0.0) then
        varm(i,k) = 0.5*(varm(i-1,k) + varm(i+1,k))

    endif


return
end subroutine terrain_filt_horiz!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine terrain_filt_vert(varm)
    use declarations
    implicit none

    real, dimension(nx,nz) :: varm
    
    if (terrain(i,k) .eq. 0.0 .and. terrain(i,k-1) .eq. 1.0) then !if the point below a physical point is inside the terrain mask... apply zero gradient condition 
        varm(i,k-1) = varm(i,k)
                
    endif

return
end subroutine terrain_filt_vert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine physical_grid_locs

use declarations

implicit none

	! Create a subroutine to assign the values/grid locations to xu, xs, zw, zs
    ! x=0.0 is the western (left) edge of the physical domain, and z=0.0 is the bottom of the physical domain

    do i = 1, nx
        xu(i) = (real(i) - 2.0) * dx ! Guarantees that when we are at u-point 2 (on our grid), we are at the western boundary (x=0.0)
        xs(i) = (real(i) - 1.5) * dx
    end do

    do k = 1, nz
        zw(k) = (real(k) - 2.0) * dz ! Gives us a lower boundary of height = 0
        zs(k) = (real(k) - 1.5) * dz
    end do

return
end subroutine physical_grid_locs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_sounding

use declarations

implicit none

! Now, open the input sounding file and read in information to fill the variables we just declared, for each level

    open(12, file='/gpfs_common/share02/lackmann/amanda/mea712finalproj/input_sounding_ivan_12z.txt')

    read(12, *, end=999) inp ! Read in surface pressure, and store it in 'inp' variable

    ! Initialize counter kk for the levels read in
    kk = 1

    ! The below line marks the start of the loop; this point is returned to at the end of the loop (goto 900), until an error or end of file
900 continue

    ! At the end of the file, it jumps to label 999, which closes the file
    ! Now read in the file and all its variables, and store them
    read(12, *, err=900, end=999) & 
        inz(kk), inth(kk), inqv(kk), inu(kk)
    inzagl(kk) = inz(kk) - inz(1) ! Gets height AGL by subtracting the first value of z from the current value of z
    inqv(kk) = inqv(kk) / 1000.0 ! Convert input sounding MR from g/kg to kg/kg 
    kk = kk + 1 ! Updates kk (counter) for the next iteration through the loop
    !write(*,*) inzagl(kk)
    goto 900 ! Repeats loop for each level until the end of the file or an error

999 close(12)
    nlevs = kk - 1

    write(*, *) 'nlevs=', nlevs
    write(*,*) 'kk=', kk

    !print *, 'inqv test', inqv(1)
    ! Now, need to interpolate values read in above to the levels where variables are actually represented in the model (grid points)

    ! First, do the near-surface points...
    ! Set the values of the variables at the lowest point (in the physical domain, z=0) to equal the values at the lowest level in the sounding

    thb(2) = inth(1) ! Setting the mean pot. temp at vertical grid point 2 (k=2) to be the first potential temp value from the sounding (kk=1)
    qb(2) = (inqv(1)) !/ 1000.0 ! Convert input sounding MR from g/kg to kg/kg
    ub(2) = inu(1)

    ! Now, interpolate the remaining grid points that fall inside the physical model boundaries
    ! Essentially, just assigning base state values to each model level from the input sounding, using interpolation
    ! Nested loop of indices for vertical grid points and indices for levels read in from sounding
    ! The if statement says that zs(k) should be in between heights kk and kk+1, so we can interpolate between heights at those values to get zs(k)

    do k=1, 2
        print *, qb(k) !testing qb
    enddo

    do k = 3, nz-1
        do kk = 1, nlevs
            if (inzagl(kk) .le. zs(k) .and. inzagl(kk+1) .ge. zs(k)) then
                weight = (inzagl(kk+1) - zs(k)) / (inzagl(kk+1) - inzagl(kk))
                thb(k) = weight * inth(kk) + (1 - weight) * inth(kk+1)
                qb(k) = (weight * inqv(kk) + (1 - weight) * inqv(kk+1)) !/ 1000.0
                print *, qb(k)
                ub(k) = weight * inu(kk) + (1 - weight) * inu(kk+1)
            endif
        enddo
    enddo

    !add initial condition for tracer field (psi)
    !if height AGL is between 1km and 2km, then initialize psi as 1, elsewhere init it as zero
    !use inzagl for height agl, inzagl is in meters, so 1000 to 2000m

    ! Check if z is between 1 km and 2 km
    !do 10000 and 2000m since our height vars are in meters
    do k = 1, nz
        do i = 1, nx
            if (zs(k) .ge. 1000.0 .and. zs(k) .le. 2000.0) then
                psi(i,k) = 1.0 !put tracer in
            else
                psi(i,k) = 0.0 !do not put tracer in
            endif
        enddo
    enddo

    open(101, file='qb_check.txt')
    write(101, *) '# k      z(m)     qb(g/kg)'

    do k = 2, nz-1
        write(101, *) '# ', k, zs(k), qb(k)
    end do

    print *, 'QVTEST,', qb(k-1)
    

return
end subroutine read_sounding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_basestate

use declarations

implicit none
    !THIS SUBROUTINE COMPUTES THE BASE STATE VARS (from input data) AND APPLIES BCs TO BASE STATE VARS


    ! Compute non-dim physical pressure at the surface (surface pressure) (pisfc)
    ! inp is in hpa, need to convert to pa (*100)

    pisfc = ((inp * 100) / 100000.0) ** (rd / cp)

    ! Now, compute some other values that we need for each vertical grid level

    do k = 2, nz-1 !can do this bc filled in BCs for thb
        thvb(k) = thb(k) * (1 + (0.61 * qb(k)))
        
    enddo
    
    ! Compute pib on the scalar levels, using the hydrostatic equation

    ! First compute for k=2, which is 0.5*dz above the physical ground
    pib(2) = pisfc - g * 0.5 * dz / (cp * thvb(2))

    ! Now for the rest of the levels, integrate upward using avg thvb for each layer
    do k = 3, nz-1
        thvbavg = 0.5 * (thvb(k) + thvb(k-1))
        pib(k) = pib(k-1) - g * dz / (cp * thvbavg) ! Integrate the hydrostatic equation as we move up through vertical levels
    enddo

    ! Compute rhos on the scalar levels

    !rhos(2) = (100000.0* (pib(2))**(cv/rd))/(rd*thvb(2))

    do k = 2, nz-1
        rhos(k) = (100000.0 * (pib(k)) ** (cv / rd)) / (rd * thvb(k))
    enddo

    ! Interpolate rhow on the w-levels

    rhow(2) = (100000.0 * pisfc ** (cv / rd)) / (rd * thvb(2))

    do k = 3, nz-1
        rhow(k) = 0.5 * (rhos(k) + rhos(k-1))
    enddo

    !BOUNDARY CONDITION FOR BASE STATE VARS

	!zero gradient condition for base state vars
	thb(1) = thb(2)
	qb(1) = qb(2)
	thb(nz) = thb(nz-1)
	qb(nz) = qb(nz-1)
	pib(nz) = pib(nz-1)
	pib(1) = pib(2) !needed for when we do the surface pressure in the skewT (doesnt interact with model, physics, but makes plot nicer)
	ub(1) = ub(2)
	ub(nz) = ub(nz-1)

    rhos(1) = rhos(2)
    rhow(1) = rhow(2)
    rhow(nz) = rhow(nz-1)
    rhos(nz) = rhos(nz-1)
    thvb(1) = thvb(2)
    thvb(nz) = thvb(nz-1)
	! Turn non-dim pressure back into physical pressure in Pa (BASE STATE)
    !p(k) is my base state physical pressure

    do k = 1, nz
        !p(k) = 100000.0 * exp((cp / rd) * log(pib(k)))
        p(k) = (pib(k)**(cp/rd))*100000.0
    enddo

    ! Model will use the total (not perturbation) winds on each grid cell, so set the initial value of u(i,k) on each level to be the same as ub(k) on that level

    do k = 1, nz-1
        do i = 1, nx-1
            u(i,k) = ub(k)
        end do
    end do

	! WE HAVE NOW SET UP THE GRID AND INITIALIZED BASE STATE VARIABLES
    ! NOW DO A QUALITY CHECK (print out the values of the base state variables we have computed)

    open(100, file='basestate_check.txt')
    write(100, *) '# k      z(m)      thb(k)      pib      pb(Pa)      rhos(kg/m^3)      qb(kg/kg)      thvb(K)      uwind(m/s)'
    
    do k = 1, nz
        write(100, *) '# ', k, zs(k), thb(k),  pib(k),  p(k),  rhos(k), qb(k), thvb(k), ub(k)
    end do
    
    close(100)

return
end subroutine compute_basestate
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeout_basestate
    use declarations
    implicit none
          
    ! Before we write out the values, need to linearly interpolate the scalar u and w points to the scalar grid point 
    ! We do this rather than leaving them at their staggered locations, so that we have all fields co-located to make plotting and calculations easier

    do i = 1, nx-1
        do k = 1, nz
            uinterp(i,k) = 0.5 * (u(i,k) + u(i+1,k))
            
            ! Need to add a boundary condition (for point nx and point nz, these formulas would overrun the array ends)
            uinterp(nx,k) = uinterp(nx-1,k)  
        end do
    end do
	
	!do sep loop for w and do nz-1

	do i=1, nx
		do k=1, nz-1
			winterp(i,k) = 0.5 * (w(i,k) + w(i,k+1))
			winterp(i,nz) = winterp(i,nz-1)
		enddo
	enddo
    
    ! Now actually write out the values to the files we just opened
    write(fmt, '(a,i4,a)') '(', nx, '(e10.4,1x))' ! Write out all of the x values on one line first

    ! Now write out all values at each Z (do for each variable)
    ! Even though base state variables are only defined in z, we can still write them out at every x-position, so all arrays are the same size

    do k = 1, nz
        write(10, fmt) (ub(k), i=1, nx)
        write(11, fmt) (thb(k), i=1, nx)
        write(12, fmt) (pib(k), i=1, nx)
        write(13, fmt) (qb(k), i=1, nx) 
		write(20,fmt) (p(k), i=1, nx)
        write(26, fmt) (terrain_mask(i,k), i=1, nx)
        
    end do

    return
end subroutine writeout_basestate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initial_rand_noise
    use declarations
    implicit none

    !add initial random noise
    !added in initial condition only
    !amplitude allows us to control the magnitude of the noise being inserted (can be set to zero)

    do k = 2, nz-1
        do i = 2, nx-1
            call random_number(rand)
            th(i,k) = th(i,k) + amplitude * (2.0*rand-1.0)
        enddo
    enddo

return
end subroutine initial_rand_noise!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine flux_predictive_eqns
    use declarations
    implicit none

    !compute thv (thetav') for all gridpoints 
    do k=1, nz
        do i=1, nx
            !first compute total thetav for all points
            thvtot(i,k) = (thb(k) + th(i,k)) * (1 + 0.61*(qb(k) + qv(i,k)))
            
            !now compute thvprime
            thv(i,k) = thvtot(i,k) - thvb(k)
            !print *, thv(i,k)
        enddo
    enddo
    !write(*,*) (th(i,k))

    !calculate vt (condensate velocities for fallout)
    do k=1, nz
        do i=1,nx
            vt(i,k) = 0.0
            if(qc(i,k) .gt. 0.001) then 
                vt(i,k) = 7.0 + 7.0 * tanh(trigpi * (qc(i,k) - 0.006)/0.020)
            endif
        enddo
    enddo

    !QC for fallout (subtract the interpolated vt value (calciulated above) from w to get ground relative flux)
    do k =2, nz-1
        do i=2, nx-1

            qcp(i,k)=qcm(i,k) - 0.5*d2t*(u(i+1,k)*(qc(i+1,k)+qc( i ,k)) & 
                        -u( i ,k)*(qc( i ,k)+qc(i-1,k)))/dx & 
                        - 0.5*d2t*( rhow(k+1)*(w(i,k+1) & 
                        -0.5*(vt(i,k+1)+vt(i,k))) & 
                        *(qc(i,k+1)+qc(i, k )) & 
                        - rhow( k )*(w(i, k ) & 
                        -0.5*(vt(i,k)+vt(i,k-1))) & 
                        *(qc(i, k )+qc(i,k-1))) &
                                                         /(dz*rhos(k))
        enddo
    enddo

    !w calculation for fallout of precip (adjust buoyancy term to simulate hydrometeor loading)
    ! do k=3, nz-1
    !     do i=2, nx-1 
    !         wp(i,k)=wm(i,k) - d2t*0.25*(((u(i+1,k)+u(i+1,k-1)) & 
    !                             *(w(i+1,k)+w(i,k)) & 
    !                             -((u(i,k)+u(i,k-1))*(w(i,k)+w(i-1,k))))/dx) & 
    !                             -0.25*d2t*(rhos(k)*((w(i,k+1)+w(i,k))**2) & 
    !                             -rhos(k-1)*((w(i,k)+w(i,k-1))**2))/(rhow(k)*dz) & 
    !                             -d2t*cp*(((thvb(k)+thvb(k-1))/(2))*((pi(i,k)-pi(i,k-1))/dz)) & 
    !                             +d2t*g*0.5*((thv(i,k)/thvb(k)+thv(i,k-1)/thvb(k-1)) - 0.5 * (qc(i,k) + qc(i,k-1)))
    !     enddo
    !     ! if (thvb(k-1) ==0.0) then 
    !     !     print *, "ERRORRERROR"
    !     ! endif

    ! enddo

    !w calculation for fallout of precip (adjust buoyancy term to simulate hydrometeor loading)
    do k=3, nz-1
        do i=2, nx-1
            wp(i,k)=wm(i,k) - 0.25*d2t*( (u(i+1,k)+u(i+1,k-1)) & 
                    *(w(i+1,k)+w( i , k )) & 
                    - (u( i ,k)+u( i ,k-1)) & 
                    *(w( i ,k)+w(i-1, k )))/dx & 
                    - 0.25*d2t*(rhos( k )*(w(i,k+1)+w(i,k))**2 & 
                    -rhos(k-1)*(w(i,k)+w(i,k-1))**2) & 
                    /(dz*rhow(k)) & 
                    - d2t*cp*0.5*(thvb(k)+thvb(k-1)) & 
                    *(pi(i,k)-pi(i,k-1))/dz & 
                    + d2t*g*( & 
                    0.5*(thv(i,k)/thvb(k)+thv(i,k-1)/thvb(k-1)) & 
                    -0.5*(qc(i,k)+qc(i,k-1))) 
        enddo
    enddo


    ! !w equation/ vert wind calculation loop
    ! do k=3, nz-1
    !     do i=2, nx-1 
    !         wp(i,k)=wm(i,k) - d2t*0.25*(((u(i+1,k)+u(i+1,k-1)) &
    !                             *(w(i+1,k)+w(i,k)) &
    !                             -((u(i,k)+u(i,k-1))*(w(i,k)+w(i-1,k))))/dx) &
    !                             -0.25*d2t*(rhos(k)*((w(i,k+1)+w(i,k))**2) &
    !                             -rhos(k-1)*((w(i,k)+w(i,k-1))**2))/(rhow(k)*dz) &
    !                             -d2t*cp*(((thvb(k)+thvb(k-1))/(2))*((pi(i,k)-pi(i,k-1))/dz)) &
    !                             +d2t*g*((thv(i,k)/thvb(k)+thv(i,k-1)/thvb(k-1))/(2))
    !     enddo
    !     ! if (thvb(k-1) ==0.0) then 
    !     !     print *, "ERRORRERROR"
    !     ! endif

    ! enddo
     
    !u equation/u-wind calculation loop
    ! do k=2,nz-1
    !     do i=2, nx-1
    !         up(i,k) = um(i,k) - d2t * 0.25 * ((u(i+1,k) + u(i,k))**2 - (u(i-1,k) + u(i,k))**2)/dx & 
    !                   - d2t * 0.25 * (rhow(k+1) * (w(i,k+1) + w(i-1,k+1)) * (u(i,k+1) + u(i,k)) - rhow(k) * (w(i,k) + w(i-1,k)) & 
    !                   * (u(i,k-1) + u(i,k)))/(rhos(k)*dz) & 
    !                   - d2t * cp * thvb(k) * (pi(i,k) - pi(i-1,k))/dx  
    !     enddo
    ! enddo

    do k=2, nz-1
        do i=2, nx-1
            up(i,k)=um(i,k) - d2t*0.25*((u(i+1,k)+u(i,k))**2 & 
                    -(u(i-1,k)+u(i,k))**2)/dx & 
                    - d2t*0.25*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1)) & 
                    *(u(i,k+1)+u(i , k )) & 
                    -rhow( k )*(w(i, k )+w(i-1, k )) & 
                    *(u(i,k-1)+u(i , k ))) & 
                    /(rhos(k)*dz) & 
                    - d2t*cp*thvb(k)*(pi(i,k)-pi(i-1,k))/dx
        enddo
    enddo

!I was having issues with theta' equation, so Cam gave me the debugging idea of splitting it up into indiv terms, which worked well
!I plan to do this for all the equations (easier to debug), but stuck with just the problem equations for now for time sake
!theta equation/theta' calc loop
    do k=2,nz-1
        do i=2,nx-1
            th_one = -(0.5/dx) * ((u(i+1,k) * (th(i+1,k) + th(i,k))) - (u(i,k) * (th(i,k) + th(i-1,k))))
        
            th_two = -(0.5/(rhos(k)*dz)) * (rhow(k+1) * w(i,k+1) * (th(i,k+1) + th(i,k)) & 
                        -(rhow(k)*w(i,k) * (th(i,k) + th(i,k-1))))

            th_three = -(0.5/dz) * (w(i,k+1) * (thb(k+1) - thb(k)) & 
                        +(w(i,k)*(thb(k) - thb(k-1))))

            !combine all terms in equation to get predicted value
            thp(i,k) = thm(i,k) + d2t * (th_one + th_two + th_three)
        enddo
    enddo

    !ADD IN SURFACE FLUXES (canned tendency to the theta' field at bottom model level)
    !mimicking a physical process, so do this right after theta' calculation
    
    if (surface_fluxes == 1) then
        do i=1,nx
            call random_number(rand)
            thp(i,2) = thp(i,2) + fluxmag * d2t * (1.0 + snoise * (2.0 * rand - 1.0))
        enddo
    endif

    !pert pressure calculation loop
    do k=2,nz-1
        do i=2,nx-1
                pip(i,k)=pim(i,k) - d2t*((cs**2)/(rhos(k)*cp*(thvb(k)**2))) &
                *(((rhos(k)*thvb(k)*(u(i+1,k)-u(i,k)))/dx) &
                +0.5*((rhow(k+1)*w(i,k+1)*(thvb(k+1)+thvb(k)) &
                -rhow(k)*w(i,k)*(thvb(k)+thvb(k-1)))/dz))
        enddo
    enddo
    
    !did same thing here as did above for theta' equation (idea from Cam) 
    !(since qv was going awry first, the LH was amplifying errors in my theta' equation)
    !qv (perturbation water vapor mixing ratio) calculation loop
    do k=2,nz-1
        do i=2,nx-1
            qv_one = -(0.5/dx) * ((u(i+1,k) * (qv(i+1,k) + qv(i,k))) - (u(i,k) * (qv(i,k) + qv(i-1,k))))
            
            qv_two = -(0.5/(rhos(k)*dz)) * (rhow(k+1) * w(i,k+1) * (qv(i,k+1) + qv(i,k)) &                  
                        -(rhow(k)*w(i,k) * (qv(i,k) + qv(i,k-1))))
    
            qv_three = -(0.5/dz) * (w(i,k+1) * (qb(k+1) - qb(k)) &                                        
                        +(w(i,k) * (qb(k) - qb(k-1))))
    
            !combine all the terms and get the prediced qv
            qvp(i,k) = qvm(i,k) + d2t * (qv_one + qv_two + qv_three)
        enddo
    enddo
    
    
    !I was running into the issue of getting condensation too early, so I split up this equation, and there were still issues,
    !then tried qv, then figured that LH would first-order amplify errors in theta' field, so broke that eqn up and it worked 
    !easier to track signs/ what you are multiply/dividing by this way
    !qc (pert total condensate mixing ratio) calculation loop
    ! do k=2,nz-1
    !     do i=2,nx-1
    !         qc_one = -(0.5/dx) * ((u(i+1,k) * (qc(i+1,k) + qc(i,k)))-(u(i,k) * (qc(i,k) + qc(i-1,k))))
            
    !         qc_two = -(0.5/(rhos(k)*dz)) * (rhow(k+1) * w(i,k+1) * (qc(i,k+1) + qc(i,k)) &                  
    !                     -(rhow(k)*w(i,k) * (qc(i,k) + qc(i,k-1))))
    
    !         !combine all the terms to get predicted qc
    !         qcp(i,k) = qcm(i,k) + d2t * (qc_one + qc_two)
    !     enddo
    ! enddo

    !psi -> perturbation value of passive tracer field calculation loop
    do k=2, nz-1
        do i=2, nx-1
            psip(i,k) = psim(i,k) - d2t*(((u(i+1,k)) * 0.5*(psi(i+1,k) + psi(i,k))) - ((u(i,k)) * 0.5*(psi(i,k) + psi(i-1,k))))/dx & 
                       - d2t * (((rhow(k+1)) * (w(i,k+1)) * 0.5*(psi(i,k+1) + psi(i,k))) - ((rhow(k)) * (w(i,k)) * 0.5*(psi(i,k) & 
                       + psi(i,k-1))))/(rhos(k)*dz)
        enddo
    enddo


    !calculates pert pressure at every grid point
    do k=1, nz
        do i=1,nx
            !ppert(i,k) = 100000.0 * (pib(k) + pi(i,k))**(cp/rd) - p(k)
            !first compute total nondim pressure
            pitot(i,k) = pi(i,k) + pib(k)
            !calculate total dim pressure from total nondim oressure
            ptot(i,k) = (pitot(i,k)**(cp/rd))*100000.0 
            !compute pert dim pressure from total dim and base dim pressure
            ppert(i,k) = ptot(i,k) - p(k) !total minus base = pert
        enddo
    enddo

    
    write(*,*) "hello"

return
end subroutine flux_predictive_eqns!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spatial_filtering
    use declarations
    implicit none

    !compute the mixing coefficients
    kmixh = cmixh * dx * dx / d2t
    kmixv = cmixv * dz * dz / d2t


    !DO ZERO GRADIENT ALONG TERRAIN FACES (if terrain is included)
    !FIST DO FOR HORIZ (so left and right terrain faces)
    !only do the future values

    do i=2, nx-1
        do k=2, nz-1
            if (terrain_mask(i,k) .ge. 1.0) then
                call terrain_filt_horiz(um)
                call terrain_filt_horiz(wm)
                call terrain_filt_horiz(thm) 
                call terrain_filt_horiz(pim)
                call terrain_filt_horiz(psim)
                call terrain_filt_horiz(qvm)
                call terrain_filt_horiz(qcm)
            endif
        enddo
    enddo

    !apply horizontal diffusion to thp
    do i = 2, nx-1
    do k = 2, nz-1
        thp(i,k) = thp(i,k) + d2t * kmixh * & 
                   (thm(i+1,k) - 2.0 * thm(i,k) + thm(i-1,k))/(dx * dx)
    enddo
    enddo

    !apply horizontal diffusion to up
    do i = 2, nx-1
        do k = 2, nz-1
            up(i,k) = up(i,k) + d2t * kmixh * & 
                       (um(i+1,k) - 2.0 * um(i,k) + um(i-1,k))/(dx * dx)
    enddo
    enddo

    !apply horiz diffusion to wp
    do i = 2, nx-1
        do k = 2, nz-1
            wp(i,k) = wp(i,k) + d2t * kmixh * & 
                       (wm(i+1,k) - 2.0 * wm(i,k) + wm(i-1,k))/(dx * dx)
    enddo
    enddo

    !apply horiz diffusion to qvp
    do i = 2, nx-1
        do k = 2, nz-1
            qvp(i,k) = qvp(i,k) + d2t * kmixh * & 
                       (qvm(i+1,k) - 2.0 * qvm(i,k) + qvm(i-1,k))/(dx * dx)
    enddo
    enddo

    !apply horiz diffusion to qcp
    do i = 2, nx-1
        do k = 2, nz-1
            qcp(i,k) = qcp(i,k) + d2t * kmixh * & 
                       (qcm(i+1,k) - 2.0 * qcm(i,k) + qcm(i-1,k))/(dx * dx)
    enddo
    enddo

    !apply horiz diffusion to psip
    do i = 2, nx-1
        do k = 2, nz-1
            psip(i,k) = psip(i,k) + d2t * kmixh * & 
                       (psim(i+1,k) - 2.0 * psim(i,k) + psim(i-1,k))/(dx * dx)
    enddo
    enddo

    !DO ZERO GRADIENT ALONG VERTICAL TERRAIN FACES (if terrain is in use)
    do i=2, nx-1
        do k=3, nz-1
            if (terrain_mask(i,k) .ge. 1.0) then
                call terrain_filt_vert(um)
                call terrain_filt_vert(wm)
                call terrain_filt_vert(thm) 
                call terrain_filt_vert(pim)
                call terrain_filt_vert(psim)
                call terrain_filt_vert(qvm)
                call terrain_filt_vert(qcm)
            endif
        enddo
    enddo

    !apply vertical diffusion to thp
    do i = 2, nx-1
    do k = 3, nz-1 
        thp(i,k) = thp(i,k) + d2t * kmixv * & 
                   (thm(i,k+1) - 2.0 * thm(i,k) + thm(i,k-1))/(dz * dz)
    enddo
    enddo

    !apply vertical diffusion to up
    !NOTE: slightly different than for other vars since we predict the total u-wind, not just perts
    !we don't want the diffusion to damp the vertical shear in base state wind

    do i = 2, nx-1
    do k = 3, nz-1
        up(i,k) = up(i,k) + d2t * kmixv * & 
                  ((um(i,k+1) - ub(k+1)) - 2.0 * (um(i,k) - ub(k)) & 
                  +(um(i,k-1) - ub(k-1)))/(dz * dz)
    enddo
    enddo
    
    !apply vertical diffusion to wp
    do i = 2, nx-1
        do k = 3, nz-1 
            wp(i,k) = wp(i,k) + d2t * kmixv * & 
                       (wm(i,k+1) - 2.0 * wm(i,k) + wm(i,k-1))/(dz * dz)
    enddo
    enddo

    !apply vertical diffusion to qvp
    do i = 2, nx-1
        do k = 3, nz-1 
            qvp(i,k) = qvp(i,k) + d2t * kmixv * & 
                       (qvm(i,k+1) - 2.0 * qvm(i,k) + qvm(i,k-1))/(dz * dz)
    enddo
    enddo

    !apply vertical diffusion to qcp
    do i = 2, nx-1
        do k = 3, nz-1 
            qcp(i,k) = qcp(i,k) + d2t * kmixv * & 
                       (qcm(i,k+1) - 2.0 * qcm(i,k) + qcm(i,k-1))/(dz * dz)
    enddo
    enddo
    
    !apply vertical diffusion to psip
    do i = 2, nx-1
        do k = 3, nz-1 
            psip(i,k) = psip(i,k) + d2t * kmixv * & 
                       (psim(i,k+1) - 2.0 * psim(i,k) + psim(i,k-1))/(dz * dz)
    enddo
    enddo

return
end subroutine spatial_filtering


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine saturation_adjust
    use declarations
    implicit none
    print *, 'entering saturation adjustment subroutine'
    !first, need to compute the total pressure, temp, and qv
    !do this using the future values that were just created, since saturation adjustment is assumed to be instantaneous
    !can just do one big loop with all of this, instead of making all of these variables new arrays
    !if were going to make them new arrays, would do qvsat(i,k) and pres_tot(i,k), etc
    !but, since they are all in the same loop, it does the calculation each time it goes through for every gridpoint
    !this way, don't need to store them as arrays
    do i=2, nx-1
        do k=2, nz-1
            pres_tot = p0 * (pip(i,k) + pib(k)) ** (cp/rd)
            temp_tot = (thp(i,k) + thb(k)) * (pip(i,k) + pib(k))
            qv_tot = qb(k) + qvp(i,k)
            
            !now compute qvsat, phi, and dqv (change in water vapor due to evap/condens)
            qvsat = (380.0/pres_tot) * exp((17.27*(temp_tot - 273.0))/(temp_tot-36.0))
            phi = qvsat * lv * 17.27 * 237.0 / (cp * (temp_tot - 36.0)**2)
            dqv = (qvsat - qvp(i,k) - qb(k)) / (1+phi)

            !now, define when we want to condense vapor vs evap condensate

            !if we have condensate or are supersaturated, evap
            if(qcp(i,k).gt.0.0 .or. qvp(i,k) + qb(k).gt.qvsat) then 
                !if supersaturated, dqv is negative ... if subsaturated, dqv is positive
                !for subsaturated, only allow evaporation up to qc, and no more
                dqv = min(dqv, qcp(i,k)) 

                qvp(i,k) = qvp(i,k) + dqv
                qcp(i,k) = qcp(i,k) - dqv
                thp(i,k) = thp(i,k) - dqv * lv / (cp*(pip(i,k) + pib(k)))
            endif
        enddo
    enddo

print *, 'ending saturation adjustment subroutine'
return
end subroutine saturation_adjust!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine modeltop_damping
    use declarations
    implicit none

    !add rayleigh damping to model top for all variables

    !FOR THETA
    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                thp(i,k) = thp(i,k) - coef*thp(i,k)
            enddo
        endif
    enddo

    !for W 
    !(*need to change zs(k) to zw(k) in the eqn, since w points are staggered abve and below u and scalar pts
   
    do k=2,nz-1
        if (zw(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zw(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                wp(i,k) = wp(i,k) - coef*wp(i,k)
            enddo
        endif
    enddo

    !for U
    !since we are predicting the total u-wind field, and we only want to damp the perts...
    !we need to subtract out the base state, so we don't damp the base state

    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                !up(i,k) = (up(i,k) - ub(k)) - coef*(up(i,k) - ub(k))
                up(i,k) = up(i,k) - coef*(up(i,k)-ub(k))
            enddo
        endif
    enddo

    !FOR QV
    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                qvp(i,k) = qvp(i,k) - coef*qvp(i,k)
            enddo
        endif
    enddo

    !FOR QC
    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                qcp(i,k) = qcp(i,k) - coef*qcp(i,k)
            enddo
        endif
    enddo

    !FOR PSI
    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                psip(i,k) = psip(i,k) - coef*psip(i,k)
            enddo
        endif
    enddo

    !for pi'
    do k=2,nz-1
        if (zs(k) .ge. zbottom) then
            coef = raydmpc * 0.5* (1.0-cos(trigpi*(zs(k)-zbottom)/(zw(nz)-zbottom)))
            do i=2,nx-1
                pip(i,k) = pip(i,k) - coef*pip(i,k)
            enddo
        endif
    enddo

return
end subroutine modeltop_damping!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timefilter
    use declarations
    implicit none

!apply asselin time filtering to all predicted variables except pi
!do on every gridpoint
    
!THETA
do i=1, nx
    do k=1, nz
        th(i,k) = th(i,k) + epsilon * (thp(i,k) - 2* th(i,k) + thm(i,k))

        !U
        u(i,k) = u(i,k) + epsilon * (up(i,k) - 2* u(i,k) + um(i,k)) 

        !W
        w(i,k) = w(i,k) + epsilon * (wp(i,k) - 2* w(i,k) + wm(i,k))

        !qv
        qv(i,k) = qv(i,k) + epsilon * (qvp(i,k) - 2* qv(i,k) + qvm(i,k))

        !qc
        qc(i,k) = qc(i,k) + epsilon * (qcp(i,k) - 2* qc(i,k) + qcm(i,k))

        !psi
        psi(i,k) = psi(i,k) + epsilon * (psip(i,k) - 2* psi(i,k) + psim(i,k))
    enddo
enddo

return
end subroutine timefilter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine handoff_values
    use declarations
    implicit none

    !this goes inside the time loop. since both handoffs need to be done at every timestep
    !only handing off present values to past needs to be done before the loop (since we need to give the scheme past info, which we are just using the ICs for that)
    !handoff present values to past
    !handoff future values to present (this helps so that when I print out th, I am printing out the predicted val at every timestep)
 
    thm = th
    th = thp 
    
    um = u
    u = up
    
    wm = w 
    w = wp 
    
    pim = pi 
    pi = pip 

    qvm = qv
    qv = qvp

    qcm = qc
    qc = qcp

    psim = psi
    psi = psip

end subroutine handoff_values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zerogradLB_pert_var_BCs
    use declarations
    implicit none

    !ZERO GRADIENT lateral boundary conditions and zero gradient top and bottom bounds
    !BOUNDARY CONDITION FOR PERTURBATION VARS

    !BOUNDARY CONDITION FOR PERTURBATION VARS
    ! IMPLEMENT A ZERO-GRADIENT CONDITION ON THE BOTTOM AND TOP BOUNDARIES (set k=1 and k=nz values for all variables)
    !only need to loop over i for top/bottom bounds (since we are setting the vertical bounds for all horiz pts)
    do i = 1, nx-1
        ! Set the top boundary (k=nz)
        
        thp(i,nz) = thp(i,nz-1)
        up(i,nz) = up(i,nz-1)
        !w(i,nz) = w(i,nz-1)
        pip(i,nz) = pip(i,nz-1)
        qvp(i,nz) = qvp(i,nz-1)
        qcp(i,nz) = qcp(i,nz-1)
        wp(i,nz) = 0
        psip(i,nz) = psip(i,nz-1)
        
        ! Set the bottom boundary (k=1)
        
        thp(i,1) = thp(i,2)
        up(i,1) = up(i,2)
        wp(i,1) = wp(i,2)
        pip(i,1) = pip(i,2)
        qvp(i,1) = qvp(i,2)
        qcp(i,1) = qcp(i,2)
        wp(i,1) = 0
        wp(i,2) = 0
        psip(i,1) = psip(i,2)
        
    end do


!also ZERO GRADIENT conditions for lateral boundaries (lateral BCs)
do k=2, nz-1
  thp(i,k) = thp(2,k)
  thp(nx,k) = thp(nx-1,k)

  qcp(i,k) = qcp(2,k)
  qcp(nx,k) = qcp(nx-1,k)

  qvp(i,k) = qvp(2,k)
  qvp(nx,k) = qvp(nx-1,k)

  psip(i,k) = psip(2,k)
  psip(nx,k) = psip(nx-1,k)

  pip(i,k) = pip(2,k)
  pip(nx,k) = pip(nx-1,k)

  wp(i,k) = wp(2,k)
  wp(nx,k) = wp(nx-1,k)

  up(2,k) = up(3,k)
  up(1,k) = up(2,k)
  up(nx,k) = up(nx-1,k)

enddo

return
end subroutine zerogradLB_pert_var_BCs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine periodicLB_pert_var_BCs
    use declarations
    implicit none
    !Periodic lateral boundary conditions and zero gradient top and bottom bounds
    !BOUNDARY CONDITION FOR PERTURBATION VARS

        !BOUNDARY CONDITION FOR PERTURBATION VARS
        ! IMPLEMENT A ZERO-GRADIENT CONDITION ON THE BOTTOM AND TOP BOUNDARIES (set k=1 and k=nz values for all variables)
    !only need to loop over i for top/bottom bounds (since we are setting the vertical bounds for all horiz pts)
    do i = 1, nx-1
        ! Set the top boundary (k=nz)
        
        thp(i,nz) = thp(i,nz-1)
        up(i,nz) = up(i,nz-1)
        !w(i,nz) = w(i,nz-1)
        pip(i,nz) = pip(i,nz-1)
        qvp(i,nz) = qvp(i,nz-1)
        qcp(i,nz) = qcp(i,nz-1)
        wp(i,nz) = 0
        psip(i,nz) = psip(i,nz-1)
        
        ! Set the bottom boundary (k=1)
        
        thp(i,1) = thp(i,2)
        up(i,1) = up(i,2)
        wp(i,1) = wp(i,2)
        pip(i,1) = pip(i,2)
        qvp(i,1) = qvp(i,2)
        qcp(i,1) = qcp(i,2)
        wp(i,1) = 0
        wp(i,2) = 0
        psip(i,1) = psip(i,2)
        
    end do


!periodic conditions for lateral boundaries (lateral BCs)
do k=1, nz
    up(nx,k) = up(2,k) !bringing the easternmost value back over to the west
    up(1,k) = up(nx-1,k) !filling in the buffer (non-physical) cells

    wp(nx,k) = wp(2,k)
    wp(1,k) = wp(nx-1,k)

    thp(nx,k) = thp(2,k)
    thp(1,k) = thp(nx-1,k)

    pip(nx,k) = pip(2,k)
    pip(1,k) = pip(nx-1,k)

    psip(nx,k) = psip(2,k)
    psip(1,k) = psip(nx-1,k)

    qcp(nx,k) = qcp(2,k)
    qcp(1,k) = qcp(nx-1,k)

    qvp(nx,k) = qvp(2,k)
    qvp(1,k) = qvp(nx-1,k)
enddo

return
end subroutine periodicLB_pert_var_BCs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine leapfrog_setup

use declarations

implicit none

    ! LEAPFROG SCHEME NEEDS PAST AND PRESENT INFO TO MAKE PREDICTION
    ! We don’t have info from the ‘minus’ time level (before the start of the model), so make a rough approximation to get values for um and thm
    ! Set um=u for every grid point and thm=th

    do i = 1, nx-1
        do k = 1, nz-1
            um(i,k) = u(i,k)
            thm(i,k) = th(i,k)
            wm(i,k) = w(i,k)
            pim(i,k) = pi(i,k)
            psim(i,k) = psi(i,k)
            qvm(i,k) = qv(i,k)
            qcm(i,k) = qc(i,k)
        end do
    end do

return
end subroutine leapfrog_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine warmbubble

use declarations

implicit none

! PUT WARM BUBBLE IN THE MODEL (initial disturbance)

    ! Now calculate r and do the piecewise function to calculate theta perturbation at the present time (th)

    do i = 2, nx-1
        do k = 2, nz-1
            r(i,k) = sqrt(((xs(i) - xc) / rx) ** 2 + ((zs(k) - zc) / rz) ** 2)
            if (r(i,k) <= 1) then
                th(i,k) = (delT * (cos(trigpi * r(i,k)) + 1)) / 2
            else
                th(i,k) = 0.0
            end if
        end do
    end do

    !write(*, *) th

return
end subroutine warmbubble
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeout_perts
    use declarations
    implicit none

    !write(fmt, '(a,i4,a)') '(', nx, '(e10.4,1x))' ! Write out all of the x values on one line first
    
    !first, calc and write out the values of u and w interpolated to scalar points
    !need to do it right before we writeout, since want to make sure the interpolated values are not used in the model integration
    !interpolate u and w winds for writing out and plotting
    do i = 1, nx-1
        do k = 1, nz
            uinterp(i,k) = 0.5 * (u(i,k) + u(i+1,k)) !changed from up to u
            
            ! Need to add a boundary condition (for point nx and point nz, these formulas would overrun the array ends)
            uinterp(nx,k) = uinterp(nx-1,k)  
        end do
    end do
	
	!!do sep loop for w and do nz-1

	do i=1, nx
		do k=1, nz-1
            winterp(i,k) = 0.5 * (w(i,k) + w(i,k+1))
			winterp(i,nz) = winterp(i,nz-1)
		enddo
	enddo

    ! quickly calcualte total theta for plotting a stratified atmosphere (base state plus pert theta)
    !first make an array of the base state values for all timesteps (so the arrays are the same size and can add)
    do i =1, nx
        do k = 1, nz
            thtot(i,k) = thb(k) + th(i,k)
        enddo
    enddo

    ! Now, Write out the variables that are perturbations on the scalar points (we still write them out on w points)
    ! That's why we do do k=1, nz for them
    do k = 1, nz
        write(14, fmt) (th(i,k), i=1, nx)
        write(15, fmt) (uinterp(i,k), i=1, nx)
        write(16, fmt) (winterp(i,k), i=1, nx)
        write(17, fmt) (ppert(i,k), i=1, nx)
        write(18, fmt) (qv(i,k), i=1, nx)
        write(19, fmt) (qc(i,k), i=1, nx)
        write(21,fmt) (thv(i,k), i=1, nx)
        write(22,fmt) (pi(i,k), i=1, nx)
        write(23,fmt) (psi(i,k), i=1, nx)
        write(24,fmt) (thtot(i,k), i=1, nx)
        !write(25,fmt) (qv_tot(i,k), i=1, nx)
    end do

    !print out maxes and mins to check
    print *, 'MAX_th=', maxval(th), '  MIN_th=', minval(th)
    print *, 'MAX_pi=', maxval(pi), '  MIN_pi=', minval(pi)
    print *, 'MAX_u=', maxval(u), '  MIN_u=', minval(u)
    print *, 'MAX_w=', maxval(w), '  MIN_w=', minval(w)
    print *, 'MAX_qv=', maxval(qv), '  MIN_qv=', minval(qv)
    print *, 'MAX_qC=', maxval(qc), '  MIN_qC=', minval(qc)
    
    !print *, 'MAX_P=', maxval(p), '  MIN_P=', minval(p)
    !print *, 'MAX_pi=', maxval(ppert), '  MIN_pi=', minval(ppert)
    print *, 'MAX_psi=', maxval(psi), '  MIN_psi=', minval(psi)
    !print *, 'MAX_psip=', maxval(psip), '  MIN_psip=', minval(psip)

    !print *, 'kk =', kk ,'   inzagl=', inz(k-1) , ' zs = ', zs(k-1)
    !print *, 'cmix=', cmixh ,'   kmix=', kmixh
    
return
end subroutine writeout_perts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine close_files
    use declarations
    implicit none

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)

    
    write(*, *) 'Done writing files'
return
end subroutine close_files