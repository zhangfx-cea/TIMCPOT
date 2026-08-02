!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftshiftcmplx(dir,x,n)
implicit none
complex::x(n)
integer::dir,n
complex::y(n)
integer::k,i,nn
if(dir.eq.-1.or.dir.eq.1)then
else
    write(*,*)"Error from fftshift: dir is wrong"
    stop
endif
if(mod(n,2).eq.0)then
    k=int(n/2.0)
else
    k=int(n/2.0)
    if(dir.eq.1)k=k+1
endif

    nn=0
    do i=k+1,n
	nn=nn+1
	y(nn)=x(i)
    enddo
    do i=1,k
	nn=nn+1
	y(nn)=x(i)
    enddo
    x=y
end subroutine fftshiftcmplx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dft(dir,n,x,xout)
implicit none
integer::dir,n
complex::x(n),xout(n),w
integer::i,j
real,parameter::pi=acos(-1.0)

if(dir.eq.1)then
elseif(dir.eq.-1)then
else
    write(*,*)"Error: dir =1 for DFT and dir =-1 for IDFT"
    stop
endif

w=cmplx(cos(-dir*2.0*pi/n),sin(-dir*2.0*pi/n))
do i=1,n
    xout(i)=0.0;
    do j=1,n
	xout(i)=xout(i)+x(j)*w**((i-1)*(j-1))
    enddo
enddo
if(dir.eq.-1)then
    do i=1,n
	xout(i)=xout(i)/(n*1.0)
    enddo
endif
end subroutine dft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft(dir,n,x)
!dir =1  for FFT
!dir =-1 for IFFT
!n must be powers of 2
implicit none
complex::x(0:n-1)
integer::dir,n
!!!!!!!!!!!!!!!!!
real,parameter::pi=acos(-1.0)
integer::nm1,nd2,m
integer::i,j,k,l
integer::le,le2,ip
complex::exch,uri,sri
integer::pad2

if(n.ne.pad2(n)) then
    write(*,*)"fft: n=",n," not a power of 2"
    stop
endif
if(dir.eq.1)then
elseif(dir.eq.-1)then
else
    write(*,*)"Error: dir =1 for FFT and dir =-1 for IFFT"
    stop
endif
nm1=n-1
nd2=n/2
m=nint(log(real(n))/log(2.0))

j=nd2
do i=1,n-2
    if(i.lt.j)then
        exch=x(i)
        x(i)=x(j)
        x(j)=exch
    endif
    k=nd2
    do while(j.ge.k)
        j=j-k
        k=k/2
    enddo
    j=j+k
enddo

do l=1,m
    le=nint(2.0**l)
    le2=le/2
    uri=(1.0,0.0)
    sri=cmplx(cos(-dir*pi/le2),sin(-dir*pi/le2))
    do j=0,le2-1
        do i=j,nm1,le
            ip=i+le2
            exch=x(ip)*uri
            x(ip)=x(i)-exch
            x(i)=x(i)+exch
        enddo
        uri=uri*sri
    enddo
enddo
if(dir.eq.-1)then
    do i=0,n-1
        x(i)=x(i)/(n*1.0)
    enddo
endif
end subroutine fft
!!!!!!!!!!!!!!!!!!!!!!!!!
integer function pad2(n)
implicit none
integer::n
pad2 = 1
do while( pad2 < n )
        pad2 = pad2 * 2
enddo
return
end
