! ----------------------------------------------------------------
! module de résolution d'un problème de laplacien par éléments finis P1
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------
module amsta01probleme

  use amsta01maillage
  use amsta01sparse

  implicit none

  type probleme
    type(maillage), pointer :: mesh
    real(kind=8), dimension(:), pointer :: uexa, g, u, f, felim
    type(matsparse) :: p_K, p_M, p_Kelim
  end type


  contains



    ! construit un probleme à partir d'un maillage
    subroutine loadFromMesh(pb,msh)
      type(probleme), intent(inout) :: pb
      type(maillage), intent(in), target :: msh
      integer :: n, nt, i
      real(kind=8) :: x, y
      ! Pointeur
      pb%mesh => msh

      n=pb%mesh%nbNodes
      nt=pb%mesh%nbTri

      allocate(pb%uexa(n),pb%g(n),pb%f(n),pb%felim(n),pb%u(n))
      pb%uexa=0.d0
      pb%g=0.d0
      pb%f=0.d0
      pb%felim=0.d0
      pb%u=0.d0
      call sparse(pb%p_K,n,n)
      call sparse(pb%p_M,n,n)
      call sparse(pb%p_Kelim,n,n)

      do i=1,n
         ! initialisation de la solution theorique, du second membre et de la condition aux limites

         ! Pour f = 0, g = x :
         ! x=pb%mesh%coords(i,1)

         ! Pour f = -2(x(x-6)+y(y-2)), g = 0 :
         x = pb%mesh%coords(i,1)
         y = pb%mesh%coords(i,2)

         !pb%uexa(i) = x*(6-x)*y*(2-y)
         !pb%f(i) = 2*(x*(6-x)+y*(2-y))
         
         pb%uexa(i) = exp((x+y)/8)
         pb%f(i) = exp((x+y)/8)/16

         ! g est la restriction de uexa sur le bord
         if (pb%mesh%refNodes(i) == 1 .OR. pb%mesh%refNodes(i) == -3) then
            pb%g(i) = pb%uexa(i)
         end if
      end do
    end subroutine loadFromMesh


    ! assemblage des matrices de rigidité et de masse, et du second membre
    subroutine assemblage(pb)
      type(probleme), intent(inout) :: pb
      real(kind=8), dimension(2) :: s1,s2,s3
      integer, dimension(3) :: s
      real(kind=8), dimension(9) :: kel, mel
      integer :: nt, i, j, k
      real(kind=8) :: x

      nt=pb%mesh%nbTri

      do i=1,nt
        s=pb%mesh%triVertices(i,1:3)
        s1=pb%mesh%coords(s(1),1:2)
        s2=pb%mesh%coords(s(2),1:2)
        s3=pb%mesh%coords(s(3),1:2)

        kel=kelem(s1,s2,s3)
        mel=melem(s1,s2,s3)

        do j=1,3
          do k=1,3
            call addtocoeff(pb%p_K,s(j),s(k),kel(3*(j-1)+k))
            call addtocoeff(pb%p_M,s(j),s(k),mel(3*(j-1)+k))
          end do
        end do
      end do

      call sort(pb%p_K)
      call sort(pb%p_M)

      pb%f=spmatvec(pb%p_M,pb%f)
    end subroutine assemblage





    ! pseudo-élimination des conditions essentielles
    !     pb : problème sur lequel appliquer la pseudo-élimination
    !     id : numéro du domaine de bord
    subroutine pelim(pb,id,id2)

      implicit none
      
      type(probleme), intent(inout) :: pb
      integer, intent(in) :: id
      integer, intent(in), optional :: id2
      
      integer, dimension(:), pointer :: indelim
      integer :: n, nn, i, ii, j, id3
      real(kind=8) :: val
      
      pb%felim=pb%f-spmatvec(pb%p_K,pb%g)
      pb%p_Kelim=pb%p_K

      n=pb%mesh%nbNodes

      if(present(id2)) then 
         nn=count(pb%mesh%refNodes == id) + count(pb%mesh%refNodes == id2)
      else
         nn=count(pb%mesh%refNodes == id)
      end if
      
      allocate(indelim(nn))

      if(present(id2)) then
         indelim=pack((/ (i, i=1,n) /), pb%mesh%refNodes == id .OR. pb%mesh%refNodes == id2)
      else
         indelim=pack((/ (i, i=1,n) /), pb%mesh%refNodes == id)
      end if
         

      do ii=1,nn
         i=indelim(ii)
         val=coeff(pb%p_K,i,i)
         pb%felim(i)=pb%g(i)*val
         do j=1,n
            if (j /= i) then
               call delcoeff(pb%p_Kelim,i,j)
            call delcoeff(pb%p_Kelim,j,i)
          end if
        end do
      end do
end subroutine pelim





    ! calcul de la matrice de rigidité élémentaire
    function kelem(s1,s2,s3) result(kel)
      real(kind=8), dimension(:), intent(in) :: s1,s2,s3
      real(kind=8), dimension(9) :: kel
      real(kind=8) :: x12,x23,x31,y12,y23,y31, a

      x12=s1(1)-s2(1)
      x23=s2(1)-s3(1)
      x31=s3(1)-s1(1)
      y12=s1(2)-s2(2)
      y23=s2(2)-s3(2)
      y31=s3(2)-s1(2)
      a=2.d0*dabs(x23*y31-x31*y23)

      kel(1)=(x23*x23+y23*y23)/a
      kel(2)=(x23*x31+y23*y31)/a
      kel(3)=(x23*x12+y23*y12)/a
      kel(4)=kel(2)
      kel(5)=(x31*x31+y31*y31)/a
      kel(6)=(x31*x12+y31*y12)/a
      kel(7)=kel(3)
      kel(8)=kel(6)
      kel(9)=(x12*x12+y12*y12)/a
    end function kelem





    ! calcul de la matrice de masse élémentaire
    function melem(s1,s2,s3) result(mel)
      real(kind=8), dimension(:), intent(in) :: s1,s2,s3
      real(kind=8), dimension(9) :: mel
      real(kind=8) :: x12,x23,x31,y12,y23,y31, a1, a2

      ! x12=s1(1)-s2(1)
      x23=s2(1)-s3(1)
      x31=s3(1)-s1(1)
      ! y12=s1(2)-s2(2)
      y23=s2(2)-s3(2)
      y31=s3(2)-s1(2)
      a1=dabs(x23*y31-x31*y23)/12.d0
      a2=a1/2.d0

      mel(1)=a1
      mel(2)=a2
      mel(3)=a2
      mel(4)=a2
      mel(5)=a1
      mel(6)=a2
      mel(7)=a2
      mel(8)=a2
      mel(9)=a1
    end function melem





    ! calcul de la solution du problème par factorisation LU
    subroutine solveLU(pb)
      type(probleme), intent(inout) :: pb
      type(matsparse) :: L, U
      call lufact(pb%p_Kelim,L,U)
      call lusolve(L,U,pb%felim, pb%u)
    end subroutine solveLU





    ! export de la solution au format vtu pour Paraview
    !     mesh : mailllage
    !     sol : vecteur solution
    !     solexa : vecteur solution exacte
    !     fname : nom du fichier de sortie (optionel)
    !             le nom doit contenir l'extension .vtu
    subroutine saveToVtu(mesh, sol, solexa, fname)
      type(maillage), intent(in) :: mesh
      real(kind=8), dimension(mesh%nbNodes), intent(in) :: sol, solexa
      character(len=*), intent(in), optional :: fname
      character(len=300) :: filename, n1, n2, tmp
      integer :: i

      filename="sol.vtu"
      if (present(fname)) then
        filename=fname
      end if

      open(unit=19, file=filename, form='formatted', status='unknown')
      write(19,*) '<VTKFile type="UnstructuredGrid" version="0.1"  byte_order="LittleEndian">'
      write(19,*) '<UnstructuredGrid>'
      n1=computeAttributeFormat("NumberOfPoints",mesh%nbNodes)
      n2=computeAttributeFormat("NumberOfCells", mesh%nbTri)
      write(19,*) '<Piece '//trim(adjustl(n1))//' '//trim(adjustl(n2))//'>'
      write(19,*) '<PointData>'
      write(19,*) '<DataArray type="Float64" Name="u" format="ascii">'
      do i=1, mesh%nbNodes
        write(19,*) sol(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Float64" Name="uexa" format="ascii">'
      do i=1, mesh%nbNodes
        write(19,*) solexa(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Float64" Name="erreur" format="ascii">'
      do i=1, mesh%nbNodes
        write(19,*) sol(i)-solexa(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</PointData>'
      write(19,*) '<Points>'
      write(19,*) '<DataArray type="Float64" Name="Nodes" NumberOfComponents="3" format="ascii">'
      do i=1, mesh%nbNodes
        write(19,*) mesh%coords(i,:)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</Points>'
      write(19,*) '<Cells>'
      write(19,*) '<DataArray type="Int32" Name="connectivity" format="ascii">'
      do i=1, mesh%nbTriTot
        write(19,*) mesh%triVertices(i,:)-1
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Int32" Name="offsets" format="ascii">'
      do i=1, mesh%nbTriTot
        write(19,*) 3*i
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="UInt8" Name="types" format="ascii">'
      do i=1, mesh%nbTriTot
        write(19,*) 5
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</Cells>'
      write(19,*) '</Piece>'
      write(19,*) '</UnstructuredGrid>'
      write(19,*) '</VTKFile>'
      close(19)
    end subroutine saveToVtu




    ! fonction qui permet de construire une chaîne de type attribut
    ! utilisée dans saveToVtu
    function computeAttributeFormat(s,i) result(n)
      character(len=*), intent(in) :: s
      integer, intent(in) :: i
      character(len=100) :: n, istr
      write(istr,*) i
      n=s//'="'//trim(adjustl(istr))//'"'
    end function computeAttributeFormat



end module amsta01probleme
