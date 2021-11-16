C   Original Author(s):
C       C.Y. Lin, E.I. Izgorodina, D. Brittain, and K. Zhang
C
C   Original script can be found at: 
C       https://rsc.anu.edu.au/~mcoote/scripts.php
C   
C   When using this codes please using the following citation: 
C       E. I. Izgorodina , C.Y. Lin, M. L. Coote. Phys. Chem. Chem. Phys. 9, 2507 (2007)
C
C   ----------------------------------------------------------------------
C
C   ConfMaker.f
C
C   Author of changes: Samuel James Higgs: s.higgs@surrey.ac.uk: 11/09/2021
C
C   Changes to the original script include:
C   
C   - GNU Fortran gcc/7.1.0-2.28 
C   - updated example usage
C   - gaussian 16 inputfile commands
C   - updated default method: b3lyp/6-31g(d) -> m062x/6-31+g(d,p)
C   - output file extension change: .com -> .gjf
C
C   ----------------------------------------------------------------------
C
C   USAGE
C       
C       $ module load gcc/7.1.0-2.28   
C       
C       $ gfortran ConfMaker.f -O3 -o ./ConfMaker -ffixed-line-length-132  
C       
C   EXE
C
C       $ ConfMaker $mol  
C
C   REQUIREMENTS   
C
C       $mol.input $mol.zmat      
C
C   DESCRIPTION
C   
C   ConfMaker mol screens all the possible combinations of specified dihedral angles
C   require mol.input and mol.zmat where mol is the input file name 
C   NOTE: mol must not have any  '.' in it's name)
C
C   ----------------------------------------------------------------------
C
C   FORMAT for mol.input
C
C       #dependent_dih         dih#         0. 120. 240. 0. 0. 0.
C                              dih#
C                              ...
C
C
C   for example:
C
C     3       dih4       0. 180.    0.   0.   0.   0.
C             dih5
C             dih6
C     2       dih7       0. 120. 240.    0.   0.   0.
C             dih8 
C     1       dih10      0.  60.  120.  180.  240.  300.
C
C   ----------------------------------------------------------------------
C
C   The first column is the number of dependent dihedral angle, ie when you rotate dih_a, if
C   dih_b and dih_c needs to rotate simultaneously, then you have to put them into the dependent list.
C   The second column is the dihedrals needs to rotate.  
C   The third-eighth column is the degrees (in real number) you are rotating. You have to put 0. to fill all the columns.

       program ConfMaker

       !implicit none

       integer*8 i,k,l,j,nrot(100000),nsumrot,ndih,ntotal,ndep(100000),nfile,m,n
       integer*8 ln,nlines
       real*8 valdih(100000),rot(100000,6)
       character*1 abc(20)
       character*80 afile(100000),dih(100000),dihinp(100000),dihdep(100000,6)
       character*80,dihinp2(100000),tmp,line(100000),bfile,bcommand,ccommand
       character*300 asystem1,asystem2,amol,azmatfile,acommand,theory,basis,solvent
       character*2 charge,multi,calc,iodine,ts,gas
       logical ex 

       data abc(1:20)/"a","b","c","d","e","f","g","h","i","j","k","l"
     $      ,"m","n","o","p","q","r","s","t"/


       call getarg(1,amol)

       theory = "M062X"
       basis = "6-31G*"
       print *," we are going to use ", trim(theory), "/",trim(basis)," for this calculation? (Y/N)"
       read *,calc
       if (calc .eq. 'N') then
           print *,"If not, What is the level of theory and basis set for this calculation? (separated by space)"
           read *,theory, basis
       endif
       print *,"Is this a TRANSITION STATE optimisation? (Y/N)"
       read *,ts
       print *,"Is this a GAS PHASE optimisation? (Y/N)"
       read *,gas
       if (gas .eq. 'N') then
           print *,"If not, What is the solvent for this PCM/UAKS calculation?"
           read *,solvent
       endif

       open(22,file='mc',status='unknown')
       print *,"What is the CHARGE and MULTIPLICITY of the molecule? (separated by space)"
       read *,charge, multi
       close(22)

       amol = trim(amol)
       write(bfile,'(a,''.zmat'')') trim(amol)
       open(30,file=trim(bfile),status='old')

       nlines = 1
       do while(.true.)
       read(30,'(a80)',end=101) line(nlines)
       nlines = nlines + 1
       enddo
 101   close(30)

       do i=1,nlines
       if (line(i) .eq. ' ') then
           j = i
           exit
       endif
      enddo

       open(40,file='zmat.part1',status='unknown')
       if (ts .eq. 'Y') then

C      to insert additional line before route section in zmat.part1 for mem option, eg %mem=1Gb use
           write(40,'(A)') "%mem=100GB"

           write(40,'(A,A,A,A,A)') "#", trim(theory), "/",trim(basis), " INT(grid=ultrafine) 
     $     OPT=(TS,calcfc,noeigentest,maxcyc=200,z-matrix)"
           if (gas .eq. 'N') then
               if (solvent .eq. 'ea') then
                   write(40,*) 'SCRF=(PCM,Solvent=ccl4,Read)'
               else
                   write(40,*) 'SCRF=(PCM,Solvent=', trim(solvent), ',Read)'
               endif
           endif
       else

C      to insert additional line before route section in zmat.part1 for mem option, eg %mem=1Gb use
           write(40,'(A)') "%mem=100GB"

           write(40,'(A,A,A,A,A)') "#", trim(theory),  "/",trim(basis)," INT(grid=ultrafine) OPT" 
           if (gas .eq. 'N') then
               if (solvent .eq. 'ea') then
                   write(40,*) 'SCRF=(PCM,Solvent=ccl4,Read)'
               else
                   write(40,*) 'SCRF=(PCM,Solvent=', trim(solvent), ',Read)'
               endif
           endif
       endif

       write(40,*) 
       write(40,*) trim(amol), '   confsearch'
       write(40,*) 
       write(40,*) trim(charge),' ', trim(multi)
       do i=1,j
           write(40,*) trim(line(i))
       enddo
       close(40)

       open(3,file='zmat',status='unknown')
       do i=j+1,nlines-1
           write(3,*) trim(line(i))
       enddo
       close(3)


       write(bfile,'(a,''.input'')') trim(amol)
       open(2,file=trim(bfile),status='old')

 
       rot=0.
       i=1
       nsumrot=0

  10   read(2,*,end=20) ndep(i),dih(i),rot(i,1:6)
       nrot(i)=2
       nsumrot=nsumrot+2
       do k=3,6
          if(rot(i,k).ne.0.) then
             nsumrot=nsumrot+1
             nrot(i)=nrot(i)+1
          endif
       enddo

      do j=2,ndep(i) 
        read(2,*) dihdep(i,j)
      enddo

       i=i+1
       goto 10

  20   ndih=i-1

       i=1

       open(3,file='zmat',status='unknown')
  30   read(3,*,end=40) dihinp(i),valdih(i)
       i=i+1
       goto 30

  40   ntotal=i-1 

!CC       write(*,*) 'ntotal',ntotal 
       close(2)
       close(3)

       
       write(bfile,'(a,''.list'')') trim(amol)
       open(10000,file=bfile,status='unknown')
        i=1
           do j=1,nrot(1)
              write(afile(j),'(''a'',i1)') j
              write(10000,'(''a'',i1)') j
              open(10,file=afile(j),status='unknown')
              do k=1,ntotal
              select case(ndep(1))
              case(1)
             if(trim(dih(i)).eq.trim(dihinp(k))) then
                    write(10,'(a,f20.5)') trim(dihinp(k)),valdih(k)+rot(i,j)
                    
                 else
                    write(10,'(a,f20.5)') trim(dihinp(k)),valdih(k)
              endif
              case(2:6)
             n=0
             if(trim(dih(i)).eq.trim(dihinp(k))) then
                write(10,'(a,f20.5)') trim(dihinp(k)),valdih(k)+rot(i,j)
                    n=n+1
             endif
             do m=2,ndep(1)
             if(trim(dihdep(i,m)).eq.trim(dihinp(k))) then
                    write(10,'(a,f20.5)') trim(dihinp(k)),valdih(k)+rot(i,j)
                    n=n+1
             endif             
             enddo
             if(n.eq.0) then
                write(10,'(a,f20.5)') trim(dihinp(k)),valdih(k)
             endif
             end select
             enddo
           enddo
         close(10)
         nsumrot=nrot(1)      

        
        nfile=nrot(1)
        valdih=0.0d0
        do i=2,ndih
           do k=2,nrot(i)
              do j=1,nsumrot
              open(10+j,file=trim(afile(j)),status='old')
!CC              dihinp=""
              nfile=nfile+1
              write(10000,'(a,a,i1)') trim(afile(j)),trim(abc(i)),k 
              write(afile(nfile),'(a,a,i1)') trim(afile(j)),trim(abc(i)),k 
!CC              write(*,*) j,trim(afile(j)),trim(afile(nfile))
              open(100000+k,file=trim(afile(nfile)),status='unknown')
!CC Making output file
               do l=1,ntotal
                  read(10+j,*) dihinp(l),valdih(l)
                  select case(ndep(i))
                  case(1)
                  if(trim(dih(i)).eq.trim(dihinp(l))) then
                    write(100000+k,'(a,f20.5)') trim(dihinp(l)),valdih(l)+rot(i,k)
                 else
                    write(100000+k,'(a,f20.5)') trim(dihinp(l)),valdih(l)
                endif
                case(2:6)
                n=0
                if(trim(dih(i)).eq.trim(dihinp(l))) then
                    write(100000+k,'(a,f20.5)') trim(dihinp(l)),valdih(l)+rot(i,k)
                    n=n+1
                endif
                do m=2,ndep(i)
                if(trim(dihdep(i,m)).eq.trim(dihinp(l))) then
                    write(100000+k,'(a,f20.5)') trim(dihinp(l)),valdih(l)+rot(i,k)
                    n=n+1
                endif 
                enddo
                if(n.eq.0) write(100000+k,'(a,f20.5)') trim(dihinp(l)),valdih(l)
               end select
               enddo         
!CC next j
           close(10+j,status='keep')
           enddo
!CC next k
           close(100000+k)
           enddo
           nsumrot=nsumrot*nrot(i)
         enddo 
         write(*,*) 'Total Number of Conformers Generated',nsumrot 
         close(10000)
      
       open(5,file='tmp',status='unknown')
       write(5,*)
       if (gas .eq. 'N') then
           if (trim(solvent) .eq. 'ea') then
               write(5,'(''Density=0.006114'')')
               write(5,'(''EPS=6.02'')')
           endif
           write(5,'(''Radii=UAKS'')')
           write(5,*)
       endif

       close(5)

         write(bfile,'(a,''.list'')') trim(amol)
         open(4,file=bfile,status='unknown')

  50     read(4,*,end=60) azmatfile
       write(acommand,'(''cat zmat.part1 '',a,2x,''tmp >'',2x,a,''.'',a,''.gjf'' )') trim(azmatfile),trim(amol),trim(azmatfile)

       write(bcommand,'(''rm '',a)') trim(azmatfile)
       
      !write(*,*) '++++ =>',acommand
       call system(acommand)
       call system(bcommand)
       
       goto 50
  60   close(4)     
 
       write(ccommand,*) 'rm tmp zmat zmat.part1 mc' 
       call system(ccommand)

       end program 
