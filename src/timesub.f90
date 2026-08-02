!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine jday2date(year,jday,mon,day)
implicit none
integer::year,jday
integer::mon,day
integer::mday(12)
integer::total
integer::i

mday(1)=31; mday(2)=28; mday(3)=31;  mday(4)=30;  mday(5)=31;  mday(6)=30;
mday(7)=31; mday(8)=31; mday(9)=30; mday(10)=31; mday(11)=30; mday(12)=31;
if(mod(year,4).eq.0)mday(2)=29;
total=0;
do i=1,12
    total=total+mday(i)
    if(jday.le.total)then
	mon=i;
	day=mday(i)-(total-jday);
	exit
    endif
enddo

end subroutine jday2date
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tadd(tin,secs,tout)
implicit none
integer::tin(6),tout(6)
real::secs
integer::t0(6)
double precision::sec0
integer::dayyear

t0=tin
sec0=real(t0(5)*1.0d0+t0(6)/1000.0d0)
sec0=sec0+secs*1.0d0
if(sec0<0.0d0)then
    do while(sec0<0.0d0)
	sec0=sec0+60.0d0;t0(4)=t0(4)-1
	if(t0(4)<0)then
	    t0(4)=t0(4)+60;t0(3)=t0(3)-1
	endif
	if(t0(3)<0)then
	    t0(3)=t0(3)+24;t0(2)=t0(2)-1
	endif
	if(t0(2)<1)then
	    t0(1)=t0(1)-1;
	    if(mod(t0(1),4).eq.0)then
		t0(2)=366
	    else
		t0(2)=365
	    endif
	endif
    enddo
elseif(sec0>=60.0d0)then
    do while(sec0>=60.0d0)
	sec0=sec0-60.0d0;t0(4)=t0(4)+1
	if(t0(4)>59)then
	    t0(4)=t0(4)-60;t0(3)=t0(3)+1
	endif
	if(t0(3)>23)then
	    t0(3)=t0(3)-24;t0(2)=t0(2)+1
	endif
	dayyear=365
	if(mod(t0(1),4).eq.0)dayyear=366
	if(t0(2)>dayyear)then
	    t0(2)=t0(2)-dayyear;t0(1)=t0(1)+1
	endif
    enddo
endif
tout(1:4)=t0(1:4)
tout(5)=int(sec0)
tout(6)=int((sec0-tout(5)*1.0d0)*1000.0d0)

end subroutine tadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function t1t2(t1in,t2in)
implicit none
integer::t1in(6),t2in(6)
integer::t1(6),t2(6)
integer::t(6)
integer::sw

integer::ty
integer::tday
double precision::tem1,tem2,t1t2dble

t1=t1in;t2=t2in
sw=0
if(t1(1).gt.t2(1))then
    t=t1;t1=t2;t2=t;sw=1
elseif(t1(1).eq.t2(1))then
    if(t1(2).gt.t2(2))then
	t=t1;t1=t2;t2=t;sw=1
    elseif(t1(2).eq.t2(2))then
	if(t1(3).gt.t2(3))then
	    t=t1;t1=t2;t2=t;sw=1
	elseif(t1(3).eq.t2(3))then
	    if(t1(4).gt.t2(4))then
		t=t1;t1=t2;t2=t;sw=1
	    elseif(t1(4).eq.t2(4))then
		if(t1(5).gt.t2(5))then
		    t=t1;t1=t2;t2=t;sw=1
		elseif(t1(5).eq.t2(5))then
		    if(t1(6).gt.t2(6))then
			t=t1;t1=t2;t2=t;sw=1
		    endif
		endif
	    endif
	endif
    endif
endif

tday=ty(t1(1),t2(1))
tday=tday-t1(2)+t2(2)

tem1=t1(3)*3600.0d0+t1(4)*60.0d0+t1(5)*1.0d0+t1(6)/1000.0d0
tem2=t2(3)*3600.0d0+t2(4)*60.0d0+t2(5)*1.0d0+t2(6)/1000.0d0
!write(*,*)tem1,tem2,real(t1(6))/1000.0,real(t2(6))/1000.0
t1t2dble=-1.0d0*tem1+tem2
t1t2=real(tday*24.0d0*3600.0d0+t1t2dble)
if(sw.eq.1)t1t2=-t1t2
end function t1t2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function ty(year1,year2)
implicit none
integer::year1,year2
integer::y1,y2,year
integer::sw,i

y1=year1;y2=year2;
sw=0
if(y1.gt.y2)then
    year=y1;y1=y2;y2=year;sw=1
endif
ty=0;
i=y1
do while(i<y2)
    if(mod(i,4).eq.0)then
	ty=ty+366;
    else
	ty=ty+365;
    endif
    i=i+1
enddo
if(sw.eq.1)ty=-ty

end function ty
