program evtlocate
implicit none
character(len=200)::outfile,argv,stapick,modeltimetable
character(len=100)::chartem
character(len=8),allocatable::kstnmf(:)
character(len=8)::kstnmtem
character(len=500)::ib
character(len=1)::temchar
character(len=100)::evtid
integer::modelsw,outsearchresultsw
real::lon1,lon2,lat1,lat2,dep1,dep2
real::lon,lat,dep
real::lons,lats,deps
integer::lonn,latn,depn
integer::n,nn,nmax,evtn
integer::i,j,k,seln

real::delta0,deltastep,dep0,depstep,evdp0,evdpstep
integer::deltanum,depnum,evdpnum
real,allocatable::time(:,:,:)
real::reipostime

real,allocatable::stla(:),stlo(:),stel(:),picktime(:),picktimerm(:)
integer,allocatable::picksw(:)
real::pickmin
integer::pickswtotal

real::delta,baz,az
real,allocatable::searchl(:),searchlrm(:)
real::searchmin

real::abnormalpicktimebase,abnormalpicktimemax
real::synthpickdiff
integer::navailable
real,allocatable::xx(:)
real::ave,rms,sd

integer::kz(6),kzevt(6)
integer::evtmon,evtday
real::evttime

real::evtlat,evtlon,evtdep
real::evtlatrms,evtlonrms,evtdeprms
real::evtlatsd, evtlonsd, evtdepsd
real::avemin,rmsmin,sdmin
real::tx(864001),yy(864001)
integer::indexn
real::pmax
integer::imin(1),imax(1)

integer::npts,nerr
real::b,e,dt

!real::delta
integer::maxx,phasen
real::tt(60),dtdd(60),dtdh(60),dddp(60),usrc(2)
character*8 phcd(60)
character*200::premodnam,ttimesdir
real::sdep

!delta=0.05
!sdep=0.0
! call depset(sdep,usrc)
! call trtm(delta,maxx,n,tt,dtdd,dtdh,dddp,phcd)
!do i=1,n
!    write(*,*)phcd(i),tt(i)
!enddo
!stop

if(iargc().lt.14)then
    write(*,'(a,$)')"Usage: evtlocate lon1 lons lon2 lat1 lats lat2 dep1 deps dep2 picktime.info abnormalpicktimemax "
    write(*,*)"1=outseachresult|0=not model.para"
    write(*,*)"picktime.info: kstnm stlat stlon stele 0/1 pick.time"
    write(*,*)"               > eventid"
    write(*,*)"model para: 0=globalmodel iasp91|ak135"
    write(*,*)"            1=usermodel   usermodelname"
    stop
endif

call getarg(1,argv);read(argv,*)lon1
call getarg(2,argv);read(argv,*)lons
call getarg(3,argv);read(argv,*)lon2
call getarg(4,argv);read(argv,*)lat1
call getarg(5,argv);read(argv,*)lats
call getarg(6,argv);read(argv,*)lat2
call getarg(7,argv);read(argv,*)dep1
call getarg(8,argv);read(argv,*)deps
call getarg(9,argv);read(argv,*)dep2
call getarg(10,stapick);
call getarg(11,argv);read(argv,*)abnormalpicktimemax
call getarg(12,argv);read(argv,*)outsearchresultsw
call getarg(13,argv);read(argv,*)modelsw
call getarg(14,modeltimetable)

if(modelsw.eq.0)then
    if(modeltimetable.eq.'iasp91' .or. modeltimetable.eq.'ak135')then
    else
	write(*,*)"model para should be iasp91 or ak135"
	stop
    endif
    maxx=60
    call getenv('HOMETTIMES',ttimesdir)
    !ttimesdir="/home/zhangfx/prog/ttime/ttimes/";
    premodnam=ttimesdir(1:len_trim(ttimesdir))//"/"//modeltimetable(1:len_trim(modeltimetable))
    call tabin(1,premodnam)
elseif(modelsw.eq.1)then
    open(unit=10,file=modeltimetable,form='unformatted',status='old')
    read(10)evdp0,evdpstep,evdpnum
    read(10)delta0,deltastep,deltanum
    read(10)dep0,depstep,depnum
    allocate(time(evdpnum,deltanum,depnum))
    do i=1,evdpnum
	read(10)((time(i,j,k),k=1,depnum),j=1,deltanum)
	read(10)
    enddo
    close(10)
endif


open(10,file=stapick,status='unknown')
n=0;nmax=0;
do 
    read(10,*,end=100)kstnmtem
    kstnmtem=adjustl(kstnmtem)
    if(kstnmtem(1:1).ne.">")then
	n=n+1
    else
	if(n.ge.nmax)nmax=n
	n=0
    endif
enddo
100 continue
n=nmax
close(10)

npts=864001
dt=0.1;
do i=1,npts
    tx(i)=0.0+real(i-1)*dt;
enddo
b=tx(1);e=tx(npts)
 call newhdr
 call setfhv('delta',  dt,  nerr)
 call setfhv('b',       b,  nerr)
 call setfhv('e',       e,  nerr)
 call setnhv('npts', npts,  nerr)
 call setihv('iftype','itime',nerr)


allocate(kstnmf(nmax),stla(nmax),stlo(nmax),stel(nmax),picksw(nmax),picktime(nmax),picktimerm(nmax))
allocate(searchl(nmax),searchlrm(nmax))
allocate(xx(nmax))
open(10,file=stapick,status='old')
nn=0;pickmin=1000000.0;
do
    read(10,'(a)',end=200)ib;
    ib=adjustl(ib)
    if(ib(1:1).ne.">")then
	nn=nn+1
	read(ib,*)kstnmf(nn),stla(nn),stlo(nn),stel(nn),picksw(nn),picktime(nn),kz(1),kz(2),kz(3),kz(4),kz(5),kz(6)
	if(picksw(nn).eq.1 .and. picktime(nn).le.pickmin)pickmin=picktime(nn)
    else
	read(ib,*)temchar,evtid
	if(outsearchresultsw.eq.1)outfile=evtid(1:len_trim(evtid))//".PICK.RM.MINVAL";
	if(outsearchresultsw.eq.1)open(30,file=outfile,status='unknown')
	pickmin=pickmin
	pickswtotal=0
	do i=1,nn
	    picktimerm(i)=picktime(i)-pickmin
	    pickswtotal=pickswtotal+picksw(i)
	    if(outsearchresultsw.eq.1)write(30,*)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),picktime(i),pickmin,picktimerm(i)
	enddo
	if(outsearchresultsw.eq.1)close(30)
	!!!!!!!!!!!!! IF no pick, then Initial and next
	if(pickswtotal.eq.0)then
	    nn=0;pickmin=1000000.0;
	    cycle
	endif

	avemin=1000000.0;rmsmin=1000000.0;sdmin=1000000.0;
	if(outsearchresultsw.eq.1)then
	    outfile=evtid(1:len_trim(evtid))//".FAST.ARS.DAT";
	    open(20,file=outfile,status='unknown')
	endif

	lon=lon1;
	do while(lon.le.lon2)
	    lat=lat1;
	    do while(lat.le.lat2)
		dep=dep1
		do while(dep.le.dep2)
			!do n=1,2
			    !if (n==1)then
				!lon=105.236;lat=28.477;dep=20.0
			    !endif
			    !if (n==2)then
				!lon=105.179;lat=28.043;dep=0.0
			    !endif
		    ! Calculate time at each station & Find the minivalue
		    searchmin=1000000.0
		    !!!!! setup evtdep firstly
		    if(modelsw.eq.0)call depset(dep,usrc)

		    do i=1,nn
			if(picksw(i).eq.1)then
			    !!!!!!!!!!!!!!! Distance
			    !call gcarc_baz_az_d(stla(i),stlo(i),lat,lon,0,1,delta,baz,az)
			    !!! First way to find distance
			    !searchl(i)=delta*6371.0;
			    !!! Second way to find distance
			    !call dist_r1_r2_gcarc(6371.0,6371.0-dep,delta,searchl(i));

			    !!!!!!!!!!!!!!! Time
			    call gcarc_baz_az_d(stla(i),stlo(i),lat,lon,0,0,delta,baz,az)
			    if(modelsw.eq.0)then
				call trtm(delta,maxx,phasen,tt,dtdd,dtdh,dddp,phcd)
				searchl(i)=tt(1)
			    elseif(modelsw.eq.1)then
				call chaxuntime(time,evdp0,evdpstep,evdpnum,delta0,deltastep,deltanum,dep0,depstep,depnum, &
						dep,delta,0.0,reipostime)
				searchl(i)=reipostime
			    endif

			    if(searchl(i).le.searchmin)then
				searchmin=searchl(i)
				seln=i
			    endif
			endif
		    enddo
		    searchmin=searchmin
!!!!!!!!!
		    pickmin=picktime(seln)
		    do i=1,nn
			picktimerm(i)=picktime(i)-pickmin
		    enddo
!!!!!!!!!
		    !!! Obtain judge-value for removing abnormal picktime
		    abnormalpicktimebase=searchmin-pickmin
!write(*,*)lon,lat,dep
		    ! Remove the minivalue & Calculate the judge-value (IF Distance and Time, divide; IF Time and Time, subtract)
		    navailable=0;
		    do i=1,nn
			!searchlrm(i)=searchl(i)-searchmin
			!write(*,*)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),searchl(i),searchmin,searchlrm(i)
			if(picksw(i).eq.1)then
			    searchlrm(i)=searchl(i)-searchmin
			    navailable=navailable+1
			    xx(navailable)=searchlrm(i)-picktimerm(i)
!write(*,*)kstnmf(i),picktime(i)-searchl(i),xx(navailable)
			    !synthpickdiff=xx(navailable)-abnormalpicktimebase
			    !if(synthpickdiff.gt.abnormalpicktimemax)then
			    if(abs(xx(navailable)).gt.abnormalpicktimemax)navailable=navailable-1
			    !navailable=navailable-1

			    !endif
			    !write(*,*)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),searchlrm(i),picktimerm(i),xx(navailable)
			endif
		    enddo
		    if(navailable.eq.0)then
			ave=0.0;rms=0.0;sd=0.0
		    else
			call avermssd(xx,navailable,ave,rms,sd)
		    endif

		    if(navailable.eq.0 .and. sd.eq.0.0)then
			ave=3.0;rms=3.0;sd=3.0
		    endif
!write(*,*)ave,rms,sd,ave-sd*3.0,ave+sd*3.0
		    if(outsearchresultsw.eq.1)write(20,*)lon,lat,dep,ave,rms,sd
		    if(sd.ne.0.0 .and. sd.le.sdmin)then
			sdmin=sd
			evtlonsd=lon; evtlatsd=lat; evtdepsd=dep;
		    endif
		    if(rms.ne.0.0 .and. rms.le.rmsmin)then
			rmsmin=rms
			evtlonrms=lon;evtlatrms=lat;evtdeprms=dep
		    endif
			!enddo
			!stop
		    dep=dep+deps
	  	enddo
		lat=lat+lats
	    enddo
	    lon=lon+lons
	enddo

	if(outsearchresultsw.eq.1)close(20)
	!write(*,*)evtid(1:len_trim(evtid)),evtlon,evtlat,evtdep,sdmin
	yy=0.0
	if(modelsw.eq.0)call depset(evtdepsd,usrc)
	do i=1,nn
	    if(picksw(i).eq.1)then
		call gcarc_baz_az_d(stla(i),stlo(i),evtlatsd,evtlonsd,0,0,delta,baz,az)
		if(modelsw.eq.0)then
		    call trtm(delta,maxx,phasen,tt,dtdd,dtdh,dddp,phcd)
		    evttime=picktime(i)-tt(1)
		elseif(modelsw.eq.1)then
		    call chaxuntime(time,evdp0,evdpstep,evdpnum,delta0,deltastep,deltanum,dep0,depstep,depnum, &
				    evtdepsd,delta,0.0,reipostime)
		    evttime=picktime(i)-reipostime
		endif
		indexn=int(evttime/0.1)+1
		yy(indexn)=yy(indexn)+10.0;
		!call tadd(kz,evttime,kzevt)
		!write(*,*)evtdep,stla(i),stlo(i),delta,tt(1),picktime(i)-tt(1),kzevt
	    endif
	enddo

	call setnhv('nzyear',kz(1),   nerr)
	call setnhv('nzjday',kz(2),   nerr)
	call setnhv('nzhour',kz(3),   nerr)
	call setnhv('nzmin', kz(4),   nerr)
	call setnhv('nzsec', kz(5),   nerr)
	call setnhv('nzmsec',kz(6),   nerr)
	call setkhv('kstnm', "SHOOT",   nerr)
	outfile=evtid(1:len_trim(evtid))//".SHOOT.TIME.SAC";
	!call wsac0(outfile,tx,yy,nerr)

	pmax=maxval(yy(1:npts))
	imax=maxloc(yy(1:npts))
	call tadd(kz,tx(imax(1)),kzevt)
	call jday2date(kzevt(1),kzevt(2),evtmon,evtday)
	!write(*,*)pmax,imax(1),tx(imax(1)),kzevt
	call gcarc_baz_az_d(evtlatsd,evtlonsd,evtlatrms,evtlonrms,0,0,delta,baz,az)
	write(*,101)evtid(1:len_trim(evtid)),kzevt(1),evtmon,evtday,kzevt(3),kzevt(4),kzevt(5),kzevt(6), &
                    evtlonsd,evtlatsd,evtdepsd,evtlonrms,evtlatrms,evtdeprms,delta,evtdepsd-evtdeprms

	nn=0;pickmin=1000000.0;
    endif
enddo
200 continue
close(10)
101 format(A,' ',I4.4,'/',I2.2,'/',I2.2,' ',I2.2,":",I2.2,":",I2.2,".",I3.3,F11.4,F11.4,F8.3,F11.4,F11.4,F8.3,F10.6,F7.2)
deallocate(xx)
deallocate(searchl,searchlrm)
deallocate(kstnmf,stla,stlo,stel,picksw,picktime)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avermssd(x,n,ave,rms,sd)
implicit none
real::x(n),ave,rms,sd
integer::n
real::temreal

integer::i
ave=0.0;rms=0.0;
do i=1,n
    ave=ave+x(i)
    rms=rms+x(i)**2
enddo
ave=ave/real(n);
rms=rms/real(n);
rms=sqrt(rms)

sd=0.0;
do i=1,n
    temreal=x(i)-ave
    sd=sd+temreal**2
enddo

sd=sd/real(n)
sd=sqrt(sd)

end subroutine avermssd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c-----------------------------------------------
	subroutine dist_r1_r2_gcarc(r1,r2,angle,dist)
!c
!c    Calculate distance by using Cosine theorem
!c    r1,r2 are in km, angle is in radian, dist is in km
!c
!c    D=r1**2+r2**2-2*r1*r2*cos(A)
!c
	implicit none
	real::r1,r2,angle,dist
	double precision::r1dble,r2dble,r12,distdble

	r1dble=dble(r1)
	r2dble=dble(r2)
	distdble=r1dble**2+r2dble**2-2.d0*r1dble*r2dble*dcos(dble(angle))
	dist=dsqrt(distdble)
	end subroutine dist_r1_r2_gcarc


!c-----------------------------------------------------------------
!c-----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!	subroutine gcarc_baz_az
!
!--------------------------------------------------------------------
!
! INPUT:
!         stalat -- Latitude for station (-90 to 90 or -pi/2 to pi/2)
!         stalon -- Longitude for station (-180 to 180 or -pi to pi)
!         evtlat -- Latitude for event (-90 to 90 or -pi/2 to pi/2)
!         evtlon -- Longitude for event (-180 to 180 or -pi to pi)
!         swin   -- switch for  input unit. 0: degree; 1: radian
!         swout  -- switch for output unit. 0: degree; 1: radian
!
! OUTPUT:
!         delta -- distance between station and event
!           baz -- azimuth for event
!            az -- azimuth for station
!
!--------------------------------------------------------------------
!
!
	subroutine gcarc_baz_az_f(stalat,stalon,evtlat,evtlon,swin,swout,delta,baz,az)
	implicit none
	real::stalat,stalon,evtlat,evtlon
	integer::swin,swout
	real::delta,baz,az
!!!!!!!!!!!!!!!!!!!!!!!
	real::stalatswp,stalonswp,evtlatswp,evtlonswp
	real,parameter::pi=acos(-1.0)
	real::piby2,rad,sph,scolat,ecolat,slon,elon
	real::a,b,c,d,e,g,h,k
	real::aa,bb,cc,dd,ee,gg,hh,kk,rhs1,rhs2

	if(stalat.eq.evtlat.and.stalon.eq.evtlon)then
	    delta=0.0;baz=0.0;az=0.0;
	    return
	endif
	piby2=pi/2.0
	rad=pi/180.0
	sph=1.0/298.257
	if(swin.eq.0)then
	    stalatswp=stalat*rad
	    stalonswp=stalon*rad
	    evtlatswp=evtlat*rad
	    evtlonswp=evtlon*rad
	elseif(swin.eq.1)then
	    stalatswp=stalat
	    stalonswp=stalon
	    evtlatswp=evtlat
	    evtlonswp=evtlon
	endif

	scolat=piby2-atan((1.0-sph)*(1.0-sph)*tan(stalatswp))
	ecolat=piby2-atan((1.0-sph)*(1.0-sph)*tan(evtlatswp))
	slon=stalonswp
	elon=evtlonswp

	a=sin(scolat)*cos(slon)
	b=sin(scolat)*sin(slon)
	c=cos(scolat)
	d=sin(slon)
	e=-cos(slon)
	g=-c*e
	h=c*d
	k=-sin(scolat)

	aa=sin(ecolat)*cos(elon)
	bb=sin(ecolat)*sin(elon)
	cc=cos(ecolat)
	dd=sin(elon)
	ee=-cos(elon)
	gg=-cc*ee
	hh=cc*dd
	kk=-sin(ecolat)

	delta=acos(a*aa+b*bb+c*cc)

	rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc-2.0
	rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k)-2.0
	baz=atan2(rhs1,rhs2)
	if(baz.lt.0.0)baz=baz+2*pi

	rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c-2.0
	rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk)-2.0
	az=atan2(rhs1,rhs2)
	if(az.lt.0.0)az=az+2*pi

	if(swout.eq.0)then
	    delta=delta/rad
	    baz=baz/rad
	    az=az/rad
	endif
	end subroutine gcarc_baz_az_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine gcarc_baz_az_d(stalat0,stalon0,evtlat0,evtlon0,swin,swout,delta0,baz0,az0)
	implicit none
	real::stalat0,stalon0,evtlat0,evtlon0
	real::delta0,baz0,az0
	integer::swin,swout
	real*8::stalat,stalon,evtlat,evtlon
	real*8::delta,baz,az
!!!!!!!!!!!!!!!!!!!!!!!
	real*8::stalatswp,stalonswp,evtlatswp,evtlonswp
	real*8,parameter::pi=dacos(-1.0d0)
	real*8::piby2,rad,sph,scolat,ecolat,slon,elon
	real*8::a,b,c,d,e,g,h,k
	real*8::aa,bb,cc,dd,ee,gg,hh,kk,rhs1,rhs2

	stalat=stalat0;stalon=stalon0;evtlat=evtlat0;evtlon=evtlon0
	if(stalat.eq.evtlat.and.stalon.eq.evtlon)then
	    delta=0.0d0;baz=0.0d0;az=0.0d0;
	    delta0=delta;baz0=baz;az0=az
	    return
	endif
	piby2=pi/2.0d0
	rad=pi/180.0d0
	sph=1.0d0/298.257d0
	if(swin.eq.0)then
	    stalatswp=stalat*rad
	    stalonswp=stalon*rad
	    evtlatswp=evtlat*rad
	    evtlonswp=evtlon*rad
	elseif(swin.eq.1)then
	    stalatswp=stalat
	    stalonswp=stalon
	    evtlatswp=evtlat
	    evtlonswp=evtlon
	endif

	scolat=piby2-datan((1.0d0-sph)*(1.0d0-sph)*dtan(stalatswp))
	ecolat=piby2-datan((1.0d0-sph)*(1.0d0-sph)*dtan(evtlatswp))
	slon=stalonswp
	elon=evtlonswp

	a=dsin(scolat)*dcos(slon)
	b=dsin(scolat)*dsin(slon)
	c=dcos(scolat)
	d=dsin(slon)
	e=-dcos(slon)
	g=-c*e
	h=c*d
	k=-dsin(scolat)

	aa=dsin(ecolat)*dcos(elon)
	bb=dsin(ecolat)*dsin(elon)
	cc=dcos(ecolat)
	dd=dsin(elon)
	ee=-dcos(elon)
	gg=-cc*ee
	hh=cc*dd
	kk=-dsin(ecolat)

	delta=dacos(a*aa+b*bb+c*cc)

	rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc-2.0d0
	rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k)-2.0d0
	baz=datan2(rhs1,rhs2)
	if(baz.lt.0.0d0)baz=baz+2.0d0*pi

	rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c-2.0d0
	rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk)-2.0d0
	az=datan2(rhs1,rhs2)
	if(az.lt.0.0d0)az=az+2.0d0*pi

	if(swout.eq.0)then
	    delta=delta/rad
	    baz=baz/rad
	    az=az/rad
	endif
	delta0=delta;baz0=baz;az0=az
	end subroutine gcarc_baz_az_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine chaxuntime(time,evdp0,evdpstep,evdpnum,delta0,deltastep,deltanum,dep0,depstep,depnum, &
              evdp,delta,reipos,reipostime)
implicit none
real::time(evdpnum,deltanum,depnum)
real::evdp0,evdpstep,delta0,deltastep,dep0,depstep,evdp,delta,reipos,reipostime
integer::evdpnum,deltanum,depnum
real::evdpactual,deltaactual,depactual
real::evdpratio,deltaratio,depratio
integer::evdpindex,deltaindex,depindex

!real::v(8)
!real r1,r2,r3,r4
!real w1,w2
if(evdp.lt.evdp0   .or.  evdp.gt.(evdp0+(evdpnum-1)*evdpstep))then
    write(*,*)"error from subroutine chaxun: evdp is out"
    stop
endif
if(delta.lt.delta0 .or.  delta.gt.(delta0+(deltanum-1)*deltastep))then
    write(*,*)"error from subroutine chaxun: delta is out"
    stop
endif
if(reipos.lt.dep0  .or.  reipos.gt.(dep0+(depnum-1)*depstep))then
    write(*,*)"error from subroutine chaxun: dep is out"
    stop
endif
evdpindex=int((evdp-evdp0)/evdpstep)+1
deltaindex=int((delta-delta0)/deltastep)+1
depindex=int((reipos-dep0)/depstep)+1

evdpactual=evdp0+(evdpindex-1)*evdpstep
deltaactual=delta0+(deltaindex-1)*deltastep
depactual=dep0+(depindex-1)*depstep

evdpratio=(evdp-evdpactual)/evdpstep
deltaratio=(delta-deltaactual)/deltastep
depratio=(reipos-depactual)/depstep
!write(*,*)evdpratio,deltaratio,depratio
 call array3d_interpolation(time,evdpnum,deltanum,depnum,evdpindex,evdpratio,deltaindex,deltaratio,depindex,depratio,reipostime)


end subroutine chaxuntime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine array3d_interpolation(array3d,nx,ny,nz,inx,rx,iny,ry,inz,rz,vout)
implicit none
real array3d(nx,ny,nz),rx,ry,rz,vout
integer nx,ny,nz,inx,iny,inz

real v1,v2,v3,v4,v5,v6,v7,v8
real r1,r2,r3,r4
real w1,w2

v1=array3d(inx  ,iny  ,inz)
v2=array3d(inx+1,iny  ,inz)
v3=array3d(inx  ,iny+1,inz)
v4=array3d(inx+1,iny+1,inz)
v5=array3d(inx  ,iny  ,inz+1)
v6=array3d(inx+1,iny  ,inz+1)
v7=array3d(inx  ,iny+1,inz+1)
v8=array3d(inx+1,iny+1,inz+1)

r1=v1*(1.0-rz)+v5*rz
r2=v2*(1.0-rz)+v6*rz
r3=v3*(1.0-rz)+v7*rz
r4=v4*(1.0-rz)+v8*rz

w1=r1*(1.0-ry)+r3*ry
w2=r2*(1.0-ry)+r4*ry
	
vout=w1*(1.0-rx)+w2*rx

end subroutine array3d_interpolation
