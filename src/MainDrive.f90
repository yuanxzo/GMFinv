!  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>
!  
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <https://www.gnu.org/licenses/>


Program maindrive
    use commondata
    Integer::info,i
    Logical::VALL1,VALL2,VALL3
    Integer(kind = 4)::time1,time2
    Character(len = 2000)::Fomt1,Fomt2,Fomt3,Fomt4

	Write(*,*) "<GMFinv>  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>"
    Write(*,*) "This program comes with ABSOLUTELY NO WARRANTY."
    Write(*,*) "This is free software, and you are welcome to redistribute it"
    Write(*,*) "under certain conditions; see <https://www.gnu.org/licenses/> for details."
    Write(*,*)
	Write(*,*)
    Write(*,*) '**********************************************************************'
    Write(*,*) '           Multimodal dispersion curves inversion program             '
    Write(*,*) '                based on generalized misfit function                  '
    Write(*,*)
    Write(*,*) 'Timestamp: 2021-03-01'
    Write(*,*) 'Version  : 1.0'
    Write(*,*) 'Brief    : '
    Write(*,*) 'Reference: '
    Write(*,*) '**********************************************************************'
    Write(*,*)
    
    
    
    !---------- Check input files ----------！
    Inquire(file='input-pars.txt',EXIST=VALL1)
    Inquire(file='input-data.txt',EXIST=VALL2)
    Inquire(file='input-mods.txt',EXIST=VALL3)
    If( .not.(VALL1 .and. vall2 .and. vall3) ) Then
        Write(*,*) "The 'input-' files in the current directory is incomplete. Press 'Enter' to end the program!"
        Write(*,*) '======================'
        Write(*,*) '* Inverted by GMFinv *'
        Write(*,*) 'E-mail: yuanxzo@qq.com'
        Write(*,*) 'Regards, Bo Yang.'
        Write(*,*) '======================='
        info=0
        Goto 999
    Else
        info=1
    End If

    !---------- Data preprocessing before inversion ----------！
    ! 1. Import particle swarm optimization parameters
    ! 2. Import the measured Rayleigh wave dispersion data, and perform point classification processing
    ! 3. Import initial model
    ! 4. Open memory for some data
    Call pretreatment

    Call Random_seed()

    Call system_clock(time1)
    
    !---------- Start inversion ----------!
    write(Fomt1,'("(","(A4),(f17.7)",I0,"(f15.7)",")")') lenchrom
    write(Fomt2,'("(",I0,"(f15.7)",")")') 4
    write(Fomt3,'("(","(I4),(I16),(E24.7)",I0,"(f13.3)",")")') lenchrom
    write(Fomt4,'("(","(I4),(f17.7)",I0,"(f15.7)",")")') lenchrom
    Write(*,'(A19)') ' In iteration......'
    Write(*,'(A41,I4,A7)') ' >>The program will perform the inversion',cyc,' times.'
    If (cyc>1) Then
        Write(*,'(A84)') ' >>In this inversion, the details of the iteration are not displayed and preserved, '
        Write(*,'(A53)') '   and only the results of each inversion are output.'
    End If

    Open(30,file='out_profile.inv',   STATUS='replace')
    Open(50,file='out_iteration.inv', STATUS='replace')
    Open(40,file='out_finalmodel.inv',STATUS='replace')
    Do icyc=1,cyc
        
        Call PSO

        If (cyc==1) Then
            Write(*,'(a81)') "The misfit value and S-wave velocity profile of inversion results is as follows:"
            Write(* ,Fomt1) '1',bestf(gen),bestone
            Write(30,Fomt1) '1',bestf(gen),bestone

            Do i=1,ceng-1
                Write(40,Fomt2) bestone(ceng-1+i),vp(i)/vs(i)*bestone(ceng-1+i),dns(i),bestone(i)
            End do
            Write(40,Fomt2) bestone(lenchrom),vp(ceng)/vs(ceng)*bestone(lenchrom),dns(ceng),0.0d0
            
            !Write(50,*) 'iter','bestin','GMF','bestone(Column 4 -> End)'
            Do i=0,gen
                Write(50,Fomt3) i,int(bestx(i,1)),bestf(i),bestx(i,2:2*ceng)
            End Do
        Else
            Write(* ,Fomt4) icyc,bestf(gen),bestone
            Write(30,Fomt4) icyc,bestf(gen),bestone
        End If
    End Do

    !---------- Total inversion time ----------!
    Call system_clock(time2)
    Print*,"The inversion took",(time2-time1)/real(1000)," s"!

    !---------- Clean cache ----------！
    Deallocate(obsf,obsv)
    Deallocate(vs,vp,dns,thk,db,ub)
    Deallocate(bestone,bestf,bestx)
    Deallocate(x,v)
	Close(30)
	Close(40)
	Close(50)
	
    !---------- End of program ----------！
    Write(*,*)
    Write(*,*) '======================'
    Write(*,*) '* Inverted by GMFinv *'
    Write(*,*) 'E-mail: yuanxzo@qq.com'
    Write(*,*) 'Regards, Bo Yang.'
    Write(*,*) '======================'

    !Read(*,*)

999 If (info==0) Then
        Read(*,*)
    End If
End Program