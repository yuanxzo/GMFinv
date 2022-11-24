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


subroutine pretreatment
    use commondata
	Use omp_lib
    Integer::i,t,k,nthreads,info
    Character::modeltype(1),nowtype(1)
    Real(8)::time1,time2
    
    !---- Read the inversion parameters
    Open(15,file='input-pars.txt')
    Read(15,*) gen
    Read(15,*) pop
    Read(15,*) bal
    Read(15,*) c1
    Read(15,*) c2
    Read(15,*) iw
    Read(15,*) GMF_Kind
    Read(15,*) cyc
    Read(15,*) nthreads
    If ( GMF_Kind>4 .or. GMF_Kind<1 ) Then
        GMF_kind=1
        Print*, 'Error: GMF_Kind must be an integer between 1 and 4! Set it to 1.'
    End If
    If ( nthreads<1 ) Then
        nthreads=1
        Print*, 'Error: nthreads must be an integer greater than 0! Set it to 1.'
    End If
    !$ call omp_set_num_threads(nthreads)
    Close(15)
    
    !---- Read measured dispersion points
    Open(20,file='input-data.txt')
    Read(20,*) nf,bo(1),bo(2)
    Allocate(obsf(nf),obsv(nf),weit(nf))
    Do i=1,nf
        Read(20,*) obsf(i),obsv(i),weit(i)
    End Do
    Close(15)
	
	!----Point classification processing
    Call PCT    

    !---- Read the initial model and set the search range of PSO algorithm
    Open(25,file='input-mods.txt')
    Read(25,*) ceng,modeltype(1)
    lenchrom=2*ceng-1
    If (bal==1) Then
        pop=10*lenchrom
    End If
    Allocate(vs(ceng),vp(ceng),dns(ceng),thk(ceng-1),db(2*ceng-1),ub(2*ceng-1))
    Do i=1,ceng
        If (i<ceng) Then
            Read(25,*) vs(i),vp(i),dns(i),thk(i)
        Else
            Read(25,*) vs(i),vp(i),dns(i)
        End If
    End Do
    Read(25,*) nowtype(1)
    Do while( nowtype(1).ne.modeltype(1) )
        Read(25,*) nowtype(1)
    End Do
    fw=0
    If (modeltype(1) == 'A') Then
        Do i=1,ceng
            If (i<ceng) Then
                Read(25,*) fw
                db(i)       =thk(i)*(1-fw(3)); ub(i)       =thk(i)*(1+fw(4));
                db(i-1+ceng)=vs(i) *(1-fw(1)); ub(i-1+ceng)=vs(i) *(1+fw(2));
            Else
                Read(25,*) fw(1),fw(2)
                db(i-1+ceng)=vs(i) *(1-fw(1)); ub(i-1+ceng)=vs(i) *(1+fw(2));
            End If
        End Do

    ElseIf (modeltype(1) == 'B') Then
        Do i=1,ceng
            If (i<ceng) Then
                Read(25,*) db(i-1+ceng),ub(i-1+ceng),db(i),ub(i)
            Else
                Read(25,*) db(i-1+ceng),ub(i-1+ceng)
            End If
        End Do
    ElseIf (modeltype(1) == 'C') Then
        Read(25,*) fw(1),fw(2)
        db(1:ceng-1)=thk(1:ceng-1)*(1-fw(1));db(ceng:2*ceng-1)=vs(1:ceng)*(1-fw(1));
        ub(1:ceng-1)=thk(1:ceng-1)*(1+fw(2));ub(ceng:2*ceng-1)=vs(1:ceng)*(1+fw(2));
    End If
    Close(25)
	
    !--- Open memory for related variables
    Allocate(bestone(2*ceng-1),bestf(0:gen),bestx(0:gen,1:2*ceng))
    Allocate(x(pop,2*ceng-1),v(pop,2*ceng-1))
    
    !--- The three empirical normalization parameters are initially set to 1
    normB=1   
    normN=1
    normM=1
End subroutine pretreatment
