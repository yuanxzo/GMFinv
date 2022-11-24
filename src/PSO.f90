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


Subroutine PSO
    Use commondata
    Integer::i,j,bin(1),t
    Real(kind=8)::rand(2),fit(pop),gbestfit,gbestone(lenchrom),pbestone(pop,lenchrom),pbestfit(pop)
    Character(len = 2000)::Fomt1

    ! Initialize the position and velocity of the particles of the PSO
    Do i=1,pop
        Do j=1,lenchrom
            Call Random_number(rand)
			rand(2)=2.0d0*rand(2)-1.0d0
            x(i,j)=rand(1)*(ub(j)-db(j))+db(j)
            v(i,j)=rand(2)*(ub(j)-db(j))*0.2d0
        End Do
    End Do
    ! Initialize GMF and find the best particles
    first_itra=1
    Call GMF(fit,x,pop)
    bin=minloc(fit)
    gbestfit=fit(bin(1))
    gbestone=x(bin(1),:)
    pbestfit=fit
    pbestone=x
    bestf(0)=fit(bin(1))
    bestx(0,1)=bin(1);bestx(0,2:2*ceng)=gbestone
    
    If (cyc==1) Then
        write(Fomt1,'("(","(I4)",I0,"(f10.2)",")")') lenchrom
        write(*,'(A60)') ' >>And in this inversion will show the details of iteration.'
    End If
	
    ! Iterative optimization
    Do t=1,gen

        Do i=1,pop
            Call Random_number(rand)
            v(i,:)=iw*v(i,:)+c1*rand(1)*(pbestone(i,:)-x(i,:))+c2*rand(2)*(gbestone-x(i,:))
            Do j=1,lenchrom
                If ( abs(v(i,j)) > (ub(j)-db(j))*0.2d0 ) Then
                    Call Random_number(rand)
                    v(i,j)=v(i,j)/abs(v(i,j))*rand(2)*(ub(j)-db(j))*0.2d0
                End If
            End Do
            x(i,:)=x(i,:)+v(i,:)
            Do j=1,lenchrom
                If ( x(i,j)<db(j) ) Then
                    Call Random_number(rand)
                    x(i,j)=db(j)+rand(2)*(ub(j)-db(j))*0.5d0
                End If
                If ( x(i,j)>ub(j) ) Then
                    Call Random_number(rand)
                    x(i,j)=ub(j)-rand(2)*(ub(j)-db(j))*0.5d0
                End If
            End Do
        End Do
        
        Call GMF(fit,x,pop)
        bin=minloc(fit)
        If (gbestfit>fit(bin(1))) Then
            gbestfit=fit(bin(1))
            gbestone=x(bin(1),:)
            bestf(t)=fit(bin(1))
        Else
            bestf(t)=bestf(t-1)
        End if
        
        Do i=1,pop
            If (pbestfit(i)>fit(i)) Then
                pbestfit(i)=fit(i)
                pbestone(i,:)=x(i,:)
            End If
        End Do
        bestx(t,1)=bin(1);bestx(t,2:2*ceng)=gbestone

        If (cyc==1) Then
            write(*,Fomt1) t, gbestone
        End If
        
    End Do
    bestone=gbestone
    
End subroutine pso
    
