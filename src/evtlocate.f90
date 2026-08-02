program evtlocate
implicit none

character(len=200)::outfile,argv,stapick,modeltimetable
character(len=200)::outfile2
character(len=100)::chartem
character(len=8),allocatable::kstnmf(:)
character(len=8)::kstnmtem
character(len=500)::ib
character(len=1)::temchar
character(len=100)::evtid,evtidtem
character(len=8)::phasetype(1000)
integer::modelsw,outsearchresultsw,bootstrapsw,pickminlim
real::lon1,lon2,lat1,lat2,dep1,dep2
real::lon,lat,dep
real::lons,lats,deps
integer::lonn,latn,depn
integer::loni,latj,depk
integer::lonnmax,latnmax,depnmax
real,allocatable::lonval(:),latval(:),depval(:)
real,allocatable::deltamatrix(:,:,:),timematrix(:,:,:,:)
integer::n,nn,nmax,evtn
integer::i,j,k,seln,kk
integer::temn,temindex
real::rand

real::delta0,deltastep,dep0,depstep,evdp0,evdpstep
integer::deltanum,depnum,evdpnum
real,allocatable::time(:,:,:)
real::reipostime

real,allocatable::stla(:),stlo(:),stel(:),picktime(:),snr(:),picktimerm(:),stelcon(:)
integer::stelsw,phtype_snr_sw,inputsw
integer,allocatable::picksw(:)
real::pickmin
integer::pickswtotal,pickswrandom,pickswreduction
integer,allocatable::pickindex(:,:),pickswcopy(:)

real::delta,baz,az,deltasd,deltarms
real::baztem,aztem
real,allocatable::searchl(:),searchlrm(:)
real::searchmin
real,allocatable::searchminmatrix(:,:,:)
integer,allocatable::selnmatrix(:,:,:)

real::abnormalpicktimebase,abnormalpicktimemax
real::synthpickdiff
integer::navailable
real,allocatable::xx(:)
real::ave,rms,sd
real,allocatable::avevol(:,:,:),rmsvol(:,:,:),sdvol(:,:,:)

integer::kz(6),kzevt(6),kt(6)
integer::evtmon,evtday
real::evttime
real,allocatable::evttimelist(:)
real::evttimemin,evttimemax,shifttime
integer::evttimeminpoi,evttimemaxpoi

real::evtlat,evtlon,evtdep,mag
real::otime,datsample
real::deltadiff,evtdepdiff
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
real::t1,t2
integer::bootstrapnum
real,allocatable::bsevtlonsd(:),bsevtlatsd(:),bsevtdepsd(:)
real,allocatable::bsevtlonrms(:),bsevtlatrms(:),bsevtdeprms(:)
real::bsavemin,bsrmsmin,bssdmin
real::herrsd,verrsd,herrrms,verrrms
real::temreal
integer::temint
real::t1t2
integer::evti,evtj,evtk

bootstrapnum=10

!delta=0.05
!sdep=0.0
! call depset(sdep,usrc)
! call trtm(delta,maxx,n,tt,dtdd,dtdh,dddp,phcd)
!do i=1,n
!    write(*,*)phcd(i),tt(i)
!enddo
!stop

if(iargc().ne.18)then
    write(*,'(a,$)')"Usage: evtlocate lon1 lons lon2 lat1 lats lat2 dep1 deps dep2 "
    write(*,'(a,$)')"picktime.info stele.con(-1,0,1) switch.for.input.type(0/1/2/3) "
    write(*,*)"max.deviation pick.min.lim(4) out.result(0=no|1=search.dat|2=pick.info) model.para bootstrap(1=yes|0=no)"
    write(*,*)       "picktime.info(sw=0): kstnm stlat stlon stele 0/1 pick.time                      kztime"
    write(*,*)       "                     > eventid"
    write(*,*)       "picktime.info(sw=1): kstnm stlat stlon stele 0/1 pick.time phasetype(P|S|U) SNR kztime"
    write(*,*)       "                     > eventid"
    write(*,'(a,$)')" picktime.info(sw=2): kstnm stlat stlon stele 0/1 pick.time                      kztime"
    write(*,*)                            "otime evttime evla evlo evdp mag evtid datsample"
    write(*,*)       "                     > eventid"
    write(*,'(a,$)')" picktime.info(sw=3): kstnm stlat stlon stele 0/1 pick.time phasetype(P|S|U) SNR kztime"
    write(*,*)                            "otime evttime evla evlo evdp mag evtid datsample"
    write(*,*)       "                     > eventid"
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
call getarg(11,argv);read(argv,*)stelsw
call getarg(12,argv);read(argv,*)inputsw
call getarg(13,argv);read(argv,*)abnormalpicktimemax
call getarg(14,argv);read(argv,*)pickminlim
call getarg(15,argv);read(argv,*)outsearchresultsw
call getarg(16,argv);read(argv,*)modelsw
call getarg(17,modeltimetable)
call getarg(18,argv);read(argv,*)bootstrapsw

lonnmax=int((lon2-lon1)/lons)+10
latnmax=int((lat2-lat1)/lats)+10
depnmax=int((dep2-dep1)/deps)+10
allocate(lonval(lonnmax),latval(latnmax),depval(depnmax))

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


allocate(kstnmf(nmax),stla(nmax),stlo(nmax),stel(nmax),picksw(nmax),picktime(nmax),snr(nmax),picktimerm(nmax),stelcon(nmax))
allocate(searchl(nmax),searchlrm(nmax))
allocate(xx(nmax))
open(10,file=stapick,status='old')
nn=0;pickmin=1000000.0;
do
    read(10,'(a)',end=200)ib;
    ib=adjustl(ib)
    if(ib(1:1).ne.">")then
	nn=nn+1
	phasetype(nn)="U";snr(nn)=0.0;mag=-12345.0;datsample=0.0;
	if(inputsw.eq.0)then
	read(ib,*)kstnmf(nn),stla(nn),stlo(nn),stel(nn),picksw(nn),picktime(nn),                      kt(1),kt(2),kt(3),kt(4),kt(5),kt(6)
	elseif(inputsw.eq.1)then
	read(ib,*)kstnmf(nn),stla(nn),stlo(nn),stel(nn),picksw(nn),picktime(nn),phasetype(nn),snr(nn),kt(1),kt(2),kt(3),kt(4),kt(5),kt(6)
	elseif(inputsw.eq.2)then
	read(ib,*)kstnmf(nn),stla(nn),stlo(nn),stel(nn),picksw(nn),picktime(nn),                      kt(1),kt(2),kt(3),kt(4),kt(5),kt(6),&
		    temreal,temint,temint,temint,temint,temint,temint,evtlat,evtlon,evtdep,mag,evtidtem,datsample
	elseif(inputsw.eq.3)then
	read(ib,*)kstnmf(nn),stla(nn),stlo(nn),stel(nn),picksw(nn),picktime(nn),phasetype(nn),snr(nn),kt(1),kt(2),kt(3),kt(4),kt(5),kt(6),&
		    temreal,temint,temint,temint,temint,temint,temint,evtlat,evtlon,evtdep,mag,evtidtem,datsample
	endif

	if(nn.eq.1)then
	    kz(1)=kt(1);kz(2)=kt(2);kz(3)=kt(3);kz(4)=kt(4);kz(5)=kt(5);kz(6)=kt(6);
	else
	    if(kz(1).ne.kt(1).or.kz(2).ne.kt(2).or.kz(3).ne.kt(3).or.kz(4).ne.kt(4).or.kz(5).ne.kt(5).or.kz(6).ne.kt(6))then
		write(*,*)"Reference time is not consistent.",kz,kt;stop
	    endif
	endif
	stelcon(nn)=stel(nn)*stelsw
	if(picksw(nn).eq.1 .and. picktime(nn).le.pickmin)pickmin=picktime(nn)
    else
	read(ib,*)temchar,evtid
	if(inputsw.eq.2 .or. inputsw.eq.3)then
	    if(evtid(1:len_trim(evtid)) .ne. evtidtem(1:len_trim(evtidtem)))then
		write(*,*)"Warning: dat maybe not consistent. Please check evtid.";stop
	    endif
	endif
	if(outsearchresultsw.eq.1)then
	    outfile=evtid(1:len_trim(evtid))//".PICK.RM.MINVAL";
	    open(30,file=outfile,status='unknown')
	endif
	pickmin=pickmin
	pickswtotal=0
	do i=1,nn
	    picktimerm(i)=picktime(i)-pickmin
	    pickswtotal=pickswtotal+picksw(i)
	    if(outsearchresultsw.eq.1)write(30,*)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),picktime(i),pickmin,picktimerm(i)
	enddo
	if(outsearchresultsw.eq.1)close(30)
	! IF no pick or less than threshold, then Initial and next
	if(pickswtotal.lt.pickminlim)then
	    nn=0;pickmin=1000000.0;
	    cycle
	endif

	avemin=1000000.0;rmsmin=1000000.0;sdmin=1000000.0;
	if(outsearchresultsw.eq.1)then
	    outfile=evtid(1:len_trim(evtid))//".FAST.ARS.DAT";
	    open(20,file=outfile,status='unknown')
	endif

	! IF read dat from TIMCPOT PICK.FINAL then adjust the search range
	if(inputsw.eq.2 .or. inputsw.eq.3)then
	    !lon1=evtlon-2.0;  lon2=evtlon+2.0;
	    !lat1=evtlat-2.0;  lat2=evtlat+2.0;
	    !dep1=evtdep-30.0; dep2=evtdep+30.0;
	endif
	lon=lon1;lonn=0
	do while(lon.le.lon2)
	    lonn=lonn+1; lonval(lonn)=lon
	    lon=lon+lons
	enddo
	lat=lat1;latn=0
	do while(lat.le.lat2)
	    latn=latn+1; latval(latn)=lat
	    lat=lat+lats
	enddo
	dep=dep1;depn=0
	do while(dep.le.dep2)
	    depn=depn+1; depval(depn)=dep
	    dep=dep+deps
	enddo


	allocate(deltamatrix(lonn,latn,nn))

	!$OMP PARALLEL private(loni,latj,i,delta,baz,az)
	!$OMP DO SCHEDULE (DYNAMIC)
	do loni=1,lonn
	    do latj=1,latn

		do i=1,nn
		    call gcarc_baz_az_d(stla(i),stlo(i),latval(latj),lonval(loni),0,0,delta,baz,az)
		    deltamatrix(loni,latj,i)=delta
		enddo

	    enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	allocate(timematrix(lonn,latn,depn,nn))
	timematrix=0.0;

	if(modelsw.eq.0)then
	!call cputime(t1)
	!!$OMP PARALLEL &
	!!$OMP & private(depk,latj,loni,i,j,usrc,delta,phasen,tt,dtdd,dtdh,dddp,phcd,seln)
	!!$OMP DO SCHEDULE (DYNAMIC)
	do depk=1,depn
	    call depset(depval(depk),usrc)
	    do latj=1,latn
		do loni=1,lonn
		    do i=1,nn
			if(picksw(i).eq.1)then
			    delta=deltamatrix(loni,latj,i)
			    call trtm(delta,maxx,phasen,tt,dtdd,dtdh,dddp,phcd)
			    if(inputsw.eq.1 .or. inputsw.eq.3)then
			    do j=1,60
				if(phcd(j)(1:1).eq.phasetype(i))then
				    seln=j;exit;
				endif
			    enddo
			    endif
			    timematrix(loni,latj,depk,i)=tt(1)
			endif
		    enddo
		enddo
	    enddo
	enddo
	!!$OMP END DO
	!!$OMP END PARALLEL
	!call cputime(t2);write(*,*)t2-t1
	endif

	allocate(searchminmatrix(lonn,latn,depn),selnmatrix(lonn,latn,depn))
	searchminmatrix=10000000.0
	!$OMP PARALLEL &
	!$OMP & private(depk,latj,loni,i,searchl,delta,reipostime)
	!$OMP DO SCHEDULE (DYNAMIC)
	do depk=1,depn
	    do latj=1,latn
		do loni=1,lonn
		    
		    do i=1,nn
			if(picksw(i).eq.1)then
			    if(modelsw.eq.0)then
				searchl(i)=timematrix(loni,latj,depk,i)
			    elseif(modelsw.eq.1)then
				delta=deltamatrix(loni,latj,i)
				call chaxuntime(time,evdp0,evdpstep,evdpnum,delta0,deltastep,deltanum,dep0,depstep,depnum, &
						depval(depk),delta,stelcon(i),reipostime)
				timematrix(loni,latj,depk,i)=reipostime
				searchl(i)=reipostime
			    endif

			    if(searchl(i).le.searchminmatrix(loni,latj,depk))then
				searchminmatrix(loni,latj,depk)=searchl(i)
				selnmatrix(loni,latj,depk)=i
			    endif

			endif
		    enddo

		enddo
	    enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL


	allocate(avevol(lonn,latn,depn),rmsvol(lonn,latn,depn),sdvol(lonn,latn,depn))
	! Generate rms & sd volume
	call generateavermssd(timematrix,searchminmatrix,selnmatrix,lonn,latn,depn,picktime,nn,picksw, &
                              abnormalpicktimemax,avevol,rmsvol,sdvol)
	if(outsearchresultsw.eq.1)then
	    do depk=1,depn
		do latj=1,latn
		    do loni=1,lonn
			write(20,*)lonval(loni),latval(latj),depval(depk),avevol(loni,latj,depk),rmsvol(loni,latj,depk),sdvol(loni,latj,depk)
		    enddo
		enddo
	    enddo
	endif

	! Find minimum value for SD & RMS
	rmsmin=1000000.0;sdmin=1000000.0;
	do i=1,lonn
	    do j=1,latn
		do k=1,depn
		    if(sdvol(i,j,k).le.sdmin)then
			sdmin=sdvol(i,j,k)
			evtlonsd=lonval(i); evtlatsd=latval(j); evtdepsd=depval(k)
			evti=i;evtj=j;evtk=k;
		    endif
		    if(rmsvol(i,j,k).le.rmsmin)then
			rmsmin=rmsvol(i,j,k)
			evtlonrms=lonval(i);evtlatrms=latval(j);evtdeprms=depval(k)
		    endif
		enddo
	    enddo
	enddo

	! Bootstrap Method
	herrsd=0.0;verrsd=0.0;herrrms=0.0;verrrms=0.0
	if(bootstrapsw.eq.1)then
	    allocate(pickindex(bootstrapnum,nn),pickswcopy(nn))
	    allocate(bsevtlonsd(bootstrapnum),bsevtlatsd(bootstrapnum),bsevtdepsd(bootstrapnum))
	    allocate(bsevtlonrms(bootstrapnum),bsevtlatrms(bootstrapnum),bsevtdeprms(bootstrapnum))
	    pickswrandom=pickswtotal*0.9
	    pickswreduction=pickswtotal-pickswrandom
	    do kk=1,bootstrapnum
		! Random data
		do i=1,nn
		    !pickindex(k,i)=picksw(i)
		    pickswcopy(i)=picksw(i)
		enddo
		
		temn=0
		do while (temn.lt.pickswreduction)
		    call RANDOM_NUMBER(rand)
		    temindex=int(rand*nn)+1
		    if(temindex.ge.1.and.temindex.le.nn.and.pickswcopy(temindex).eq.1)then
			temn=temn+1
			pickswcopy(temindex)=0
		    endif
		enddo

		! Generate rms & sd volume
		!call cputime(t1)
		call generateavermssd(timematrix,searchminmatrix,selnmatrix,lonn,latn,depn,picktime,nn,pickswcopy, &
                            abnormalpicktimemax,avevol,rmsvol,sdvol)
		!call cputime(t2)
		!write(*,*)t2-t1

		bsrmsmin=1000000.0;bssdmin=1000000.0;
		do i=1,lonn
		    do j=1,latn
			do k=1,depn
			    if(sdvol(i,j,k).le.bssdmin)then
				bssdmin=sdvol(i,j,k)
				bsevtlonsd(kk)=lonval(i); bsevtlatsd(kk)=latval(j); bsevtdepsd(kk)=depval(k)
			    endif
			    if(rmsvol(i,j,k).le.bsrmsmin)then
				bsrmsmin=rmsvol(i,j,k)
				bsevtlonrms(kk)=lonval(i);bsevtlatrms(kk)=latval(j);bsevtdeprms(kk)=depval(k)
			    endif
			enddo
		    enddo
		enddo
		!write(*,*)bsevtlonsd(kk),bsevtlatsd(kk),bsevtdepsd(kk),bsevtlonrms(kk),bsevtlatrms(kk),bsevtdeprms(kk)

	    enddo
	    call generatelocationuncertainty(bsevtlonsd, bsevtlatsd, bsevtdepsd, bootstrapnum,herrsd, verrsd)
	    call generatelocationuncertainty(bsevtlonrms,bsevtlatrms,bsevtdeprms,bootstrapnum,herrrms,verrrms)

	    deallocate(pickindex,pickswcopy)
	    deallocate(bsevtlonsd,bsevtlatsd,bsevtdepsd)
	    deallocate(bsevtlonrms,bsevtlatrms,bsevtdeprms)
	endif
	! End Bootstrap Method

	deallocate(deltamatrix)
	deallocate(searchminmatrix,selnmatrix)
	deallocate(avevol,rmsvol,sdvol)

	if(outsearchresultsw.eq.1)close(20)

	allocate(evttimelist(nn))
	evttimemin=1.0e10;evttimemax=-1.0e10;
	if(modelsw.eq.0)call depset(evtdepsd,usrc)
	do i=1,nn
	    if(picksw(i).eq.1)then
		call gcarc_baz_az_d(stla(i),stlo(i),evtlatsd,evtlonsd,0,0,delta,baz,az)
		if(modelsw.eq.0)then
		    call trtm(delta,maxx,phasen,tt,dtdd,dtdh,dddp,phcd)
		    evttime=picktime(i)-tt(1)
		    evttimelist(i)=evttime
		elseif(modelsw.eq.1)then
		    call chaxuntime(time,evdp0,evdpstep,evdpnum,delta0,deltastep,deltanum,dep0,depstep,depnum, &
				    evtdepsd,delta,stelcon(i),reipostime)
		    evttime=picktime(i)-reipostime
		    evttimelist(i)=evttime
		endif
		if(evttimelist(i).le.evttimemin)evttimemin=evttimelist(i)
		if(evttimelist(i).ge.evttimemax)evttimemax=evttimelist(i)

		!indexn=int(evttime/0.1)+1
		!yy(indexn)=yy(indexn)+10.0;
		!call tadd(kz,evttime,kzevt)
		!write(*,*)evtdep,stla(i),stlo(i),delta,tt(1),picktime(i)-tt(1),kzevt
	    endif
	enddo

	! shift evttime if necessary
	if((evttimemax-evttimemin).gt.86400.0)then
	    write(*,*)"Unexpected error for ",evtid(1:len_trim(evtid))
	    nn=0;pickmin=1000000.0;
	    cycle
	endif

	if(evttimemin.ge.0.0 .and. evttimemax.le.86400.0)then
	    evttimeminpoi=0;evttimemaxpoi=0;
	else
	    evttimeminpoi=0;evttimemaxpoi=0;
	    if(evttimemin.lt.0.0)then
		evttimeminpoi=-(int(abs(evttimemin)/86400.0)+1)
	    endif
	    if(evttimemin.ge.86400.0)then
		evttimeminpoi=int((evttimemin-86400.0)/86400.0)+1
	    endif
	    if(evttimemax.le.0.0)then
		evttimemaxpoi=-(int(abs(evttimemax)/86400.0)+1)
	    endif
	    if(evttimemax.gt.86400.0)then
		evttimemaxpoi=int((evttimemax-86400.0)/86400.0)+1
	    endif
	    if(evttimeminpoi.eq.0 .and. evttimemaxpoi.eq.0)then
		write(*,*)"Error1 for ",evtid(1:len_trim(evtid));stop
	    endif
	    if(evttimeminpoi.gt.evttimemaxpoi)then
		write(*,*)"Error2 for ",evtid(1:len_trim(evtid));stop
	    endif
	endif
	shifttime=real(evttimeminpoi+evttimemaxpoi)/2.0*86400.0;
	if((evttimemin-shifttime).ge.0.0 .and. (evttimemax-shifttime).le.86400.0)then
	else
	    write(*,*)"Error3 for ",evtid(1:len_trim(evtid));stop
	endif

	! Find the shoot time
	yy=0.0
	do i=1,nn
	    if(picksw(i).eq.1)then
		evttime=evttimelist(i)-shifttime
		indexn=int(evttime/0.1)+1
		if(indexn.lt.1 .or. indexn.gt.864001)then
		    write(*,*)"Error4 for ",evtid(1:len_trim(evtid));stop
		endif
		yy(indexn)=yy(indexn)+10.0;
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
	call tadd(kz,tx(imax(1))+shifttime,kzevt)
	call jday2date(kzevt(1),kzevt(2),evtmon,evtday)
	!write(*,*)pmax,imax(1),tx(imax(1)),kzevt
	call gcarc_baz_az_d(evtlatsd,evtlonsd,evtlatrms,evtlonrms,0,0,delta,baz,az)
	write(*,101)evtid(1:len_trim(evtid)),kzevt(1),evtmon,evtday,kzevt(3),kzevt(4),kzevt(5),kzevt(6), &
                    evtlonsd,evtlatsd,evtdepsd,evtlonrms,evtlatrms,evtdeprms,delta,evtdepsd-evtdeprms,sdmin,rmsmin, &
                    herrsd,verrsd,herrrms,verrrms
	if(outsearchresultsw.eq.2)then
	    outfile=evtid(1:len_trim(evtid))//".PICK.INFO";
	    open(20,file=outfile,status='unknown')
	    outfile2=evtid(1:len_trim(evtid))//".PICK.FINAL.RELOCATION";
	    open(40,file=outfile2,status='unknown')
	    do i=1,nn
		call gcarc_baz_az_d(stla(i),stlo(i),evtlatrms,evtlonrms,0,0,deltarms,baz,az)
		call gcarc_baz_az_d(stla(i),stlo(i),evtlatsd, evtlonsd, 0,0,deltasd, baz,az)
		write(20,102)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),picktime(i),kz(1),kz(2),kz(3),kz(4),kz(5),kz(6), &
			     picktime(i)-(tx(imax(1))+shifttime), deltasd,deltarms
		otime=t1t2(kz,kzevt)
		deltadiff=0.0;evtdepdiff=0.0;
		if(inputsw.eq.2 .or. inputsw.eq.3)then
		    call gcarc_baz_az_d(evtlat,evtlon,evtlatsd,evtlonsd,0,0,deltadiff,baztem,aztem)
		    evtdepdiff=evtdepsd-evtdep
		endif
		write(40,103)kstnmf(i),stla(i),stlo(i),stel(i),picksw(i),picktime(i),phasetype(i),snr(i),kz(1),kz(2),kz(3),kz(4),kz(5),kz(6), &
			     otime, kzevt(1),kzevt(2),kzevt(3),kzevt(4),kzevt(5),kzevt(6),evtlatsd,evtlonsd,evtdepsd,mag, &
                             evtid(1:len_trim(evtid)),datsample,deltasd,baz,az,deltadiff,evtdepdiff,timematrix(evti,evtj,evtk,i)
	    enddo
	    write(40,'(A," ",A)')">",evtid(1:len_trim(evtid))
	    close(20)
	    close(40)
	endif
	deallocate(evttimelist)
	deallocate(timematrix)

	nn=0;pickmin=1000000.0;
    endif
enddo
200 continue
close(10)
101 format(A,' ',I4.4,'/',I2.2,'/',I2.2,' ',I2.2,":",I2.2,":",I2.2,".",I3.3, &
           F11.4,F11.4,F8.3,F11.4,F11.4,F8.3,F10.6,F7.2,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)
102 format(A,' ',F11.4,' ',F11.4,' ',F8.3,' ',I1.1,' ',F11.3,' ',I4.4,' ',I3.3,' ',I2.2,' ',I2.2,' ',I2.2,' ',I3.3,' ', &
           F11.3,' ',F11.4,' ',F11.4)
103 format(A,' ',F11.4,' ',F11.4,' ',F8.4,' ',I1.1,' ',F11.3,' ',A1,' ', &
           F11.3,' ',I4.4,' ',I3.3,' ',I2.2,' ',I2.2,' ',I2.2,' ',I3.3,' ', &
           F11.3,' ',I4.4,' ',I3.3,' ',I2.2,' ',I2.2,' ',I2.2,' ',I3.3,' ',F11.4,' ',F11.4,' ',F8.3,' ',F8.1,' ', &
           A,' ',F7.3,' ',F11.4,' ',F11.4,' ',F11.4,' ',F10.4,' ',F8.3,' ',F11.4)
deallocate(xx)
deallocate(searchl,searchlrm)
deallocate(kstnmf,stla,stlo,stel,stelcon,picksw,picktime,snr)
deallocate(lonval,latval,depval)
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine generatelocationuncertainty(lon,lat,dep,n,herr,verr)
implicit none
real::lon(n),lat(n),dep(n),herr,verr
integer::n

real::lonave,latave,depave
real::temreal

integer::i
real::hxx(n),vxx(n)

 call avermssd(lon,n,lonave,temreal,temreal)
 call avermssd(lat,n,latave,temreal,temreal)
 call avermssd(dep,n,depave,temreal,temreal)

do i=1,n
    call gcarc_baz_az_d(latave,lonave,lat(i),lon(i),0,0,hxx(i),temreal,temreal)
    vxx(i)=abs(depave-dep(i))
enddo
 call avermssd(hxx,n,temreal,herr,temreal)
 call avermssd(vxx,n,temreal,verr,temreal)
end subroutine generatelocationuncertainty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine generateavermssd(timematrix,searchminmatrix,selnmatrix,lonn,latn,depn,picktime,nn,pickswcopy, &
                            limitedmax,avevol,rmsvol,sdvol)
implicit none
real::timematrix(lonn,latn,depn,nn),searchminmatrix(lonn,latn,depn)
integer::selnmatrix(lonn,latn,depn)
integer::lonn,latn,depn,nn
real::picktime(nn)
integer::pickswcopy(nn)
real::limitedmax
real::avevol(lonn,latn,depn),rmsvol(lonn,latn,depn),sdvol(lonn,latn,depn)

integer::depk,latj,loni,i
real::searchmin,searchlrm(nn)
integer::seln
real::pickmin,picktimerm(nn)
integer::navailable
real::xx(nn)
real::ave,rms,sd

!$OMP PARALLEL &
!$OMP & private(depk,latj,loni,i,searchmin) &
!$OMP & private(seln,pickmin,navailable,picktimerm,searchlrm,xx,ave,rms,sd)
!$OMP DO SCHEDULE (DYNAMIC)
do depk=1,depn
    do latj=1,latn
	do loni=1,lonn
	    searchmin=searchminmatrix(loni,latj,depk)
	    seln=selnmatrix(loni,latj,depk)
	    pickmin=picktime(seln)

	    navailable=0
	    do i=1,nn
		if(pickswcopy(i).eq.1)then
		    searchlrm(i)=timematrix(loni,latj,depk,i)-searchmin
		    picktimerm(i)=picktime(i)-pickmin
		    navailable=navailable+1
		    xx(navailable)=searchlrm(i)-picktimerm(i)
		    if(abs(xx(navailable)).gt.limitedmax)navailable=navailable-1
		endif
	    enddo

	    if(navailable.eq.0)then
		ave=3.0;rms=3.0;sd=3.0
	    else
		call avermssd(xx,navailable,ave,rms,sd)
	    endif
	    avevol(loni,latj,depk)=ave
	    rmsvol(loni,latj,depk)=rms
	    sdvol(loni,latj,depk)=sd
	enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine generateavermssd
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cputime(t2)                                               
implicit none
real::t2
integer::count_1, count_rate, count_max

 call system_clock(count_1, count_rate, count_max)

t2 = (1.d0*count_1)/(1.d0*count_rate) 

return
end
