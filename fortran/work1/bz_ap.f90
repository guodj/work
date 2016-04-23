!将原始的一小时一次的数据改为一天一次。
    program dfhskhfgks
     implicit none

     integer u1,u2,i,j,year,day,h,l,apt
     character fname1*50,fname2*50
     real bze,bzm,ap,bzet,bzmt
     data u1,u2,fname1,fname2,bze,bzm,ap,l            &
          /1,2,'/home/gdj/study/data/bz_ap.txt',      &
          '/home/gdj/study/data/daybz_ap.txt',3*0.0,0/
     open(1,file=fname1,status='old')
     open(2,file=fname2,status='new')
     do while(.true.)
        do i=1,24,1
         read(1,"(i5,i4,i3,2f6.1,i4)",end=10)year,day,h,bzet,bzmt,apt
           if(bzet<=100)then
             bze=bze+bzet;bzm=bzm+bzmt
             l=l+1
           endif
           ap=ap+apt/24.0
        enddo
        if(l/=0)then
           bze=bze/l;bzm=bzm/l
        endif
        write(2,"(2i4,2f10.4,f5.1)")year,day,bze,bzm,ap
        bze=0.0;bzm=0.0;ap=0.0;l=0
     enddo
10   end program
