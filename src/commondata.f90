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
