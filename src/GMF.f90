! MIT License
! 
! Copyright (c) 2022 Bo YANG (seism.yang@foxmail.com)
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.


! This program provides four kinds of generalized objective functions
Subroutine GMF(fit,xx,nx)
    Use commondata
    Use Dispersion
    
    Integer::i,j,nx
    Real(kind=8)::Vr(Bnn,nx),Fxx(Nnn),Nfit(nx),Bfit(nx),NNfit(nx),BBfit(nx),fit(nx),xx(nx,lenchrom),dm(nx)
    Type(Model)::mods(nx)

    Do i=1,nx
        mods(i)%ceng=ceng
        mods(i)%vs(1:ceng)=xx(i,ceng:lenchrom)
        mods(i)%vp(1:ceng)=vp/vs*xx(i,ceng:lenchrom)
        mods(i)%dns(1:ceng)=dns
        mods(i)%thk(1:ceng-1)=xx(i,1:ceng-1)
    End Do

    If (GMF_Kind==1) Then  !***************************************************************GMF
        If (Bnn==0) Then
			!---------- Modified Thomoson-Haskell algorithm ----------！
			Do i=1,nx
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)       
			End Do
			fit=sqrt(Nfit)
		ElseIf (Nnn==0) Then
			Do i=1,nx	
				!---------- Fast vector-transfer algorithm ----------！
			   Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
			   Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
			End Do
			fit=sqrt(Bfit)
        Else
			Do i=1,nx
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)

				!---------- Fast vector-transfer algorithm ----------！
				Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
				Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
			End Do
			fit=sqrt(Bfit+Nfit)
        End If
        RETURN
    Elseif (GMF_Kind==2) Then  !***************************************************************GMF+dm
		If (Bnn==0) Then
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))                                          
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)      
			End Do
			fit=sqrt(Nfit+dm)
		ElseIf (Nnn==0) Then
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))                                          
				!---------- Fast vector-transfer algorithm ----------！
			    Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
			    Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
			End Do
			fit=sqrt(Bfit+dm)
		Else
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))                                          
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)
				!---------- Fast vector-transfer algorithm ----------！
				Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
				Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
			End Do
			fit=sqrt(Bfit+Nfit+dm)
        End If
        RETURN
	Elseif (GMF_Kind==3) Then  !***************************************************************NGMF
		If (Bnn==0) Then
			Do i=1,nx	
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)
				
                if (first_itra==1) Then
                    NNfit(i)=sum(Fxx**2)/real(Nnn)    
                End if
			End Do
			
			!Get the value of the empirical normalization parameter
			If (first_itra==1) Then                                                    
			    first_itra=0                                                
			    normN=abs(sum(NNfit)/nx)                                             
			End If  
			
			fit=sqrt(Nfit/normN)
		ElseIf (Nnn==0) Then
			Do i=1,nx
				!---------- Fast vector-transfer algorithm ----------！
			    Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
                Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
                if (first_itra==1) Then
                    BBfit(i)=sum((Vr(:,i)-B_v)**2)/real(Bnn)
                End if
			End Do
			
			!Get the value of the empirical normalization parameter
			If (first_itra==1) Then                                                    
			    first_itra=0                                                            
			    normB=abs(sum(BBfit)/nx)                                                
			End If 
			
			fit=sqrt(Bfit/normB)
		Else
			Do i=1,nx
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)
                
				!---------- Fast vector-transfer algorithm ----------！
				Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
				Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)

				if (first_itra==1) Then
					NNfit(i)=sum(Fxx**2)/real(Nnn)
                    BBfit(i)=sum((Vr(:,i)-B_v)**2)/real(Bnn)
                End if
			End Do
			
			!Get the values of the empirical normalization parameters
			If (first_itra==1) Then
				first_itra=0
				normB=abs(sum(BBfit)/nx)
				normN=abs(sum(NNfit)/nx)
			End If
			
			fit=sqrt(Bfit/normB+Nfit/normN)
        End If
        RETURN
    Elseif (GMF_Kind==4) Then  !***************************************************************NGMF+dm
		If (Bnn==0) Then
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))                                          
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn) 

				if (first_itra==1) Then
					NNfit(i)=sum(Fxx**2)/real(Nnn)
                End if				
			End Do
			
			!Get the values of the empirical normalization parameters
			If (first_itra==1) Then                                                     
			    first_itra=0                                                            
			    normM=abs(sum(dm)/nx)                                                   
			    normN=abs(sum(NNfit)/nx)                                                
			End If    
			
			fit=sqrt(Nfit/normN+dm/normM)
		ElseIf (Nnn==0) Then
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))                                          
				!---------- Fast vector-transfer algorithm ----------！
			    Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
				Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
				
				if (first_itra==1) Then
                    BBfit(i)=sum((Vr(:,i)-B_v)**2)/real(Bnn)
                End if
			End Do
			
			!Get the values of the empirical normalization parameters
			If (first_itra==1) Then                                                     
			    first_itra=0                                                            
			    normM=abs(sum(dm)/nx)                                                   
			    normB=abs(sum(BBfit)/nx)                                                
			End If 
			
			fit=sqrt(Bfit/normB+dm/normM)
		Else
			Do i=1,nx
				!---------- Calculate the smoothness of the model ----------！
				dm(i)=0
				Do j=1,ceng-1                                                           
				    dm(i)=dm(i)+((mods(i)%vs(j+1)-mods(i)%vs(j))/mods(i)%thk(j))**2     
				End do                                                                  
				dm(i)=sqrt(dm(i)/real(ceng-1))     
				
				!---------- Modified Thomoson-Haskell algorithm ----------！
				Do j=1,Nnn
					Call Haskell_s(Fxx(j),mods(i),N_f(j),N_v(j))
				End Do
				Nfit(i)=sum((N_w*Fxx)**2)/real(Nnn)
				
				!---------- Fast vector-transfer algorithm ----------！
				Call FVTA_c(Vr(:,i),mods(i),B_f,Bnn,1)
				Bfit(i)=sum((B_w*(Vr(:,i)-B_v))**2)/real(Bnn)
				
				if (first_itra==1) Then
					NNfit(i)=sum(Fxx**2)/real(Nnn)
                    BBfit(i)=sum((Vr(:,i)-B_v)**2)/real(Bnn)
                End if
			End Do
			
			!Get the values of the empirical normalization parameters
			If (first_itra==1) Then
				first_itra=0
				normM=abs(sum(dm)/nx)                                                   
				normB=abs(sum(BBfit)/nx)
				normN=abs(sum(NNfit)/nx)
			End If
			
			fit=sqrt(Bfit/normB+Nfit/normN+dm/normM)
        End If
        RETURN
    End IF
End subroutine GMF




! Point classification technology
Subroutine PCT
    Use commondata
    Integer::nbo,t,k
    Integer,dimension(:),allocatable::boo,Bn,Nn

    ! Set B is empty: the second value of the first line of "input-data.txt" is set to 0
    If (bo(1)==0) Then
        Bnn=0
        Nnn=nf
        Allocate(N_f(nf),N_v(nf),N_w(nf))
        N_f=obsf
        N_v=obsv
        N_w=weit
        return
    End If
    
    Bnn=bo(2)-bo(1)+1
    Nnn=nf-Bnn
    ! Set N is empty
    if (Nnn==0) Then
        Allocate(B_f(nf),B_v(nf),B_w(nf))
        B_f=obsf
        B_v=obsv
        B_w=weit
        return
    End If    

    ! Set L contains both Set B and Set N
    Allocate(B_f(Bnn),B_v(Bnn),B_w(Bnn),Bn(Bnn),N_f(Nnn),N_v(Nnn),N_w(Nnn),Nn(Nnn))
    Do i=1,Bnn
        Bn(i)=bo(1)+i-1
    End Do
    t=1;k=1
    Do i=1,nf
        If (i==Bn(t)) Then
            t=t+1
            If (t>Bnn) Then
                t=Bnn
            End If
        Else
            Nn(k)=i
            k=k+1
        End If
    End Do
    B_f=obsf(Bn)
    B_v=obsv(Bn)
    B_w=weit(Bn)
    N_f=obsf(Nn)
    N_v=obsv(Nn)
    N_w=weit(Nn)

    Deallocate(Bn,Nn)
End Subroutine PCT
