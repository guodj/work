
 program jdfsk
    implicit none
    character fname1*100,fname2*100,str(27)
    integer u1,u2,year,bday,bart,p,l,polarity(27)
    data u1,fname1,u2,fname2 &
         /1,'/home/gdj/study/graduation_project/data/initial_data/IMF_polarity_table.txt', &
          2,'/home/gdj/study/graduation_project/data/processed_data/IMF_polarity_table_nu.txt'/
    open(1,file=fname1)
    open(2,file=fname2)
    read(1,*);read(1,*)
    do while (.true.) 
        read(1,*,end=10)year,bday,bart,(str(l),l=1,27)
        do p=1,27
                selectcase(str(p))
                case('B')
                        polarity(p)=1
                case('R')
                        polarity(p)=-1
                case('W')
                        polarity(p)=9
                case('Y')
                        polarity(p)=0
                endselect
        enddo
        write(2,"(2i4,27i3)")year,bday,polarity
    enddo
10  close(u1);close(u2)
    end
