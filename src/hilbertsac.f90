program hilbertsac
implicit none
character(len=100)::argv
character(len=300)::infile,outfile,outfilenew,outdir
real::b,delta
real,allocatable::yf(:),yout(:)
complex,allocatable::youtcmplx(:)
integer::npts,nptsout,pad2,nerr
integer::nptspad
integer::i,j,n,filenum
integer::ii
integer::len0,lentrim
real::fttype
real,parameter::pi=acos(-1.0)

real,allocatable::envelope(:),phase(:),instfreq(:)
real::addval,phasediff,doublepi

doublepi=2.0*pi

if(iargc().lt.3)then
    write(*,*)"Usage:hilbertsac infiles outdir DFT(=1.0)/FFT(>1.0)"
    write(*,*)"      Notes for outfilename, e.g. if give A.sac"
    write(*,*)"      then output three files in outdir: A.sac.envelope A.sac.phase A.sac.instfreq"
    stop
endif
filenum=iargc()-2
 !call getarg(1,argv); read(argv,'(a)')infile
 call getarg(iargc()-1,argv); read(argv,'(a)')outdir
 call getarg(iargc(),argv); read(argv,*)fttype
if(fttype.lt.1.0)then
    write(*,*)"Error: DFT(=1.0)/FFT(>1.0)?"
    stop
endif

!write(*,*)outdir(1:len_trim(outdir)),fttype,filenum
do ii=1,filenum
    call getarg(ii,argv); read(argv,'(a)')infile

    open(10,file=infile, form='unformatted',access='direct',recl=4,status='old')
    read(10,rec=80)npts
    close(10)

    allocate(yf(npts),yout(npts),youtcmplx(npts))
    yf=0.0;
    call rsac1(infile,yf,npts,b,delta,npts,nerr)
    if(nerr.ne.0)then
	write(*,*)"read error";stop
    endif

    call ht(yf,yout,npts,fttype)
    !write(*,*)yout
    allocate(envelope(npts),phase(npts),instfreq(npts))

    do i=1,npts
	envelope(i)=abs(cmplx(yf(i),yout(i)))
	phase(i)=atan2(yout(i),yf(i))
	phasediff=0.0
	if(i>1)phasediff=abs(phase(i)-phase(i-1))
	if(phasediff>=pi)then
	    ! 1st step: set the diffval lt 2*pi
	    n=int(phasediff/doublepi)+1
	    addval=doublepi*real(n)
	    if(abs(phase(i)-addval-phase(i-1))<=doublepi)then
		phase(i)=phase(i)-addval
	    endif
	    if(abs(phase(i)+addval-phase(i-1))<=doublepi)then
		phase(i)=phase(i)+addval
	    endif
	    ! 2nd step: if diffval is still gt pi, then add/substract just one 2*pi
	    phasediff=abs(phase(i)-phase(i-1))
	    if(phasediff>=pi)then
		!write(*,*)i,phase(i)*180.0/pi,phase(i-1)*180.0/pi,phasediff*180.0/pi
		if(abs(phase(i)-doublepi-phase(i-1))<=pi)then
	  	    !if(i==3954)write(*,*)"OKKK1"
		    phase(i)=phase(i)-doublepi
		endif
		if(abs(phase(i)+doublepi-phase(i-1))<=pi)then
	  	    !if(i==3954)write(*,*)"OKKK2"
		    phase(i)=phase(i)+doublepi
		endif
	    endif
	    phasediff=abs(phase(i)-phase(i-1))
	    ! sometime phasediff may have little error, so great than pi
	    !if(phasediff>=pi)then
		!write(*,*)i,phase(i)*180.0/pi,phase(i-1)*180.0/pi,phasediff*180.0/pi,addval*180.0/pi;stop
	    !endif
	endif
    enddo

    do i=2,npts
	instfreq(i)=(phase(i)-phase(i-1))/delta/(2.0*pi)
    enddo
    instfreq(1)=instfreq(2)

    lentrim=len_trim(infile)
    !write(*,*)infile(1:lentrim)
    len0=1;
    do j=lentrim,1,-1
	if(infile(j:j) .eq. "/")then
	    len0=j+1; exit;
	endif
    enddo
    outfile=infile(len0:lentrim)
    !write(*,*)outfile(1:len_trim(outfile))

    b=b;
    delta=delta
    outfilenew=outdir(1:len_trim(outdir))//"/"//outfile(1:len_trim(outfile))//".envelope"
    call wsac1(outfilenew,envelope,npts,b,delta,nerr)

    outfilenew=outdir(1:len_trim(outdir))//"/"//outfile(1:len_trim(outfile))//".phase"
    call wsac1(outfilenew,phase,npts,b,delta,nerr)

    outfilenew=outdir(1:len_trim(outdir))//"/"//outfile(1:len_trim(outfile))//".instfreq"
    call wsac1(outfilenew,instfreq,npts,b,delta,nerr)

    deallocate(yf,yout,youtcmplx)
    deallocate(envelope,phase,instfreq)

enddo

end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine htcmplx(yin,yout,npts,fttype)
implicit none
real::yin(npts),fttype
complex::yout(npts)
integer::npts
complex,allocatable::yt(:),yk(:)
integer,allocatable::h(:)
integer::i,nptspad
integer::pad2

if(fttype.eq.1.0)then
    nptspad=npts
elseif(fttype.gt.1.0)then
    nptspad=pad2(int(fttype*npts))
endif
write(*,*)nptspad
allocate(yt(nptspad),yk(nptspad))
yt=0.0;
do i=1,npts
    yt(i)=cmplx(yin(i),0.0)
enddo

if(nptspad.ne.pad2(nptspad)) then
    call dft(1,nptspad,yt,yk)
else
    yk=yt
    call fft(1,nptspad,yk)
endif

allocate(h(nptspad))
do i=1,nptspad
    if(i==1 .or. i==(nptspad/2+1))then
	h(i)=1;
    elseif(i<=nptspad/2)then
	h(i)=2;
    else
	h(i)=0;
    endif
    yk(i)=yk(i)*cmplx(h(i),0.0)
enddo

if(nptspad.ne.pad2(nptspad)) then
    call dft(-1,nptspad,yk,yt)
else
    yt=yk
    call fft(-1,nptspad,yt)
endif

 !call fft(-1,nptspad,yk)
do i=1,npts
    yout(i)=yt(i)
enddo

deallocate(yt,yk,h)
end subroutine htcmplx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ht(yin,yout,npts,fttype)
implicit none
real::yin(npts),yout(npts)
integer::npts
real::fttype
complex,allocatable::yt(:),yk(:)
integer::i
integer::nptspad,nptspaddiv
integer::pad2

if(fttype.eq.1.0)then
    nptspad=npts
elseif(fttype.gt.1.0)then
    nptspad=pad2(int(fttype*npts))
endif

allocate(yt(nptspad),yk(nptspad))
yt=0.0;
do i=1,npts
    yt(i)=cmplx(yin(i),0.0)
enddo

if(nptspad.ne.pad2(nptspad)) then
    call dft(1,nptspad,yt,yk)
else
    yk=yt
    call fft(1,nptspad,yk)
endif

 call fftshiftcmplx(1,yk,nptspad)
nptspaddiv=int(nptspad/2)

do i=1,nptspaddiv
    yk(i)=yk(i)*cmplx(0.0,1.0)
enddo
do i=nptspaddiv+1,nptspad
    yk(i)=yk(i)*cmplx(0.0,-1.0)
enddo

 call fftshiftcmplx(-1,yk,nptspad)

yk(1)=0.0
if(mod(nptspad,2)==0) yk(nptspad/2+1)=0.0

if(nptspad.ne.pad2(nptspad)) then
    call dft(-1,nptspad,yk,yt)
else
    yt=yk
    call fft(-1,nptspad,yt)
endif

do i=1,npts
    yout(i)=REAL(yt(i))
enddo
deallocate(yt,yk)
end subroutine ht
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
