! MIT License
! 
! Copyright (c) 2022 Bo YANG (b.yang@petalmail.com)
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


Module commondata
    implicit none
    
    !---- The parameter that used to specify what kind of GMF is used for inversion
    Integer::GMF_Kind
    
    !---- The number of repeated inversion is cyc, and the index of the current inversion is icyc
    Integer::cyc,icyc,first_itra
     
    !---- Parameters of PSO algorithm
    Integer::bal,gen,pop
    Real(kind=8)::c1,c2,iw
    
    !---- Measured data
	!     The number of dispersion points nf, frequency vector obsf, phase velocity vector obsv, 
	!     the start and end index of the basic dispersion point bo(1) and bo(2)
	!     The weight of dispersion points weit
    !     Bnn=length(B_f .or. B_v .or. B_w),Nnn=length(N_f .or. N_v .or. N_w)
    Integer::nf,bo(2),Bnn,Nnn
    Real(kind=8),dimension(:),allocatable::obsf,obsv,weit
    Real(kind=8),dimension(:),allocatable::B_f,B_v,B_w,N_f,N_v,N_w
    
    !---- initial model
	!     The number of layers ceng, the number of inverted variables lenchrom
	!     Search range related parameters: fw, db, ub
	!     Initial model parameters: vs, vp, dns, thk
    Integer::ceng,lenchrom
    Real(kind=8)::fw(4)
    Real(kind=8),dimension(:),allocatable::vs,vp,dns,thk,db,ub
    
    !---- Inversion results and iteration curves
    Real(kind=8),dimension(:),allocatable::bestone,bestf
    Real(kind=8),dimension(:,:),allocatable::x,v,bestx
    
    !---- The three empirical normalization parameters
    Real(kind=8)::normB,normN,normM
End module commondata
