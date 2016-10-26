! ----------------------------------------------------------------
! module d'algebre lineaire sur les matrices sparse
!
! On stocke uniquement les indice de ligne et de colonne,
! ainsi que la valeur d'un coefficient non nul d'une matrice (i, j, val)
! dans 3 tableaux separes
! On stocke egalement les dimensions de la matrice, ainsi qu'un
! indicateur d'allocation des tableaux
!
! Les coefficients ne sont pas necessairement ranges dans l'ordre
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------
module amsta01sparse
  implicit none
  type matsparse
    integer :: n, m                            ! nombre de lignes et nombre de colonnes de la matrice
    integer, dimension(:), pointer :: i, j     ! tableaux des indices de ligne et de colonne des coefficients non nuls
    real(kind=8), dimension(:), pointer :: val ! tableau des valeurs des coefficients
    logical :: isAllocated                     ! true si les 3 tableaux sont alloues (a la meme taille)
  end type
  interface find
    module procedure findpos, findind
  end interface
  interface operator (+)
    module procedure spadd
  end interface
  interface operator (-)
      module procedure spminus
    end interface
    interface assignment (=)
    module procedure spaffect, spcopy
  end interface
  interface operator (*)
    module procedure spscalmat,spmatscal,spmatvec,spmatmat
  end interface
  contains



    ! initialisation d'une matrice
    !     a : matrice
    !     n : nombre de lignes
    !     m : nombre de colonnes
    !     nnz : nombre de coefficients non nuls (optionel)
    !     si nnz est precise, les tableaux sont alloues
    subroutine sparse(a, n, m, nnz)
      implicit none
      type(matsparse), intent(out) :: a
      integer, intent(in) :: n, m
      integer, intent(in), optional :: nnz
      a%n=n
      a%m=m
      a%isAllocated=.false.
      if (present(nnz)) then
        allocate(a%i(nnz),a%j(nnz),a%val(nnz))
        a%isAllocated=.true.
      end if
    end subroutine sparse




    ! allocation d'une matrice non allouee
    subroutine msallocate(a,nnz)
      type(matsparse), intent(inout) :: a
      integer, intent(in) :: nnz
      if (a%isAllocated.eqv..true.) then
        stop "ERROR in msallocate: matsparse already allocated"
      end if
      allocate(a%i(nnz),a%j(nnz),a%val(nnz))
      a%isAllocated=.true.
    end subroutine msallocate




    ! desaloccation d'une matrice
    subroutine msdeallocate(a)
      type(matsparse), intent(inout) :: a

      if (a%isAllocated.eqv..false.) then
        stop "ERROR in msdeallocate: matsparse already deallocated"
      end if
      deallocate(a%i,a%j,a%val)
      a%isAllocated=.false.
    end subroutine msdeallocate




    ! reallocation d'une matrice
    subroutine msreallocate(a,nnz)
      type(matsparse), intent(inout) :: a
      integer, intent(in) :: nnz
      if (a%isAllocated.eqv..false.) then
        stop "ERROR in msreallocate: matsparse already deallocated"
      end if
      deallocate(a%i,a%j,a%val)
      allocate(a%i(nnz),a%j(nnz),a%val(nnz))
      a%isAllocated=.true.
    end subroutine msreallocate




    ! fonction de modification/insertion d'un coefficient
    !     si la matrice est vide, elle contient cet unique coefficient
    !     sinon, soit on ecrase le coefficient s'il existe, soit on le cree
    subroutine setcoeff(a, i, j, val)
      implicit none
      integer, intent(in) :: i, j
      type(matsparse), intent(inout) :: a
      real(kind=8), intent(in) :: val
      integer, dimension(:), pointer :: itmp
      integer, dimension(:), pointer :: jtmp
      real(kind=8), dimension(:), pointer :: vtmp
      integer :: nnz, ind
      logical :: isEmpty

      if (i>a%n) then
        stop "ERROR in setcoeff: row index out of bounds" 
      end if
      if (j>a%m) then
        stop "ERROR in setcoeff: column index out of bounds"
      end if

      isEmpty=.false.
      if (a%isAllocated.eqv..false.) then
        call msallocate(a,1)
        isEmpty=.true.
      end if

      if (size(a%i) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays i and j in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%i) /= size(a%val)) then
        stop "ERROR in setcoeff: arrays i and val in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%val) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays j and val in matsparse must have the same size (did you use msallocate ?)"
      end if

      if (isEmpty.eqv..true.) then
        a%i(1)=i
        a%j(1)=j
        a%val(1)=val
      else
        ind=findpos(a, i, j)
        if (ind /= 0) then ! coeff already exists
          a%val(ind)=val
        else ! coeff does not exists, we add it to the storage
          nnz=size(a%i)
          allocate(itmp(nnz), jtmp(nnz), vtmp(nnz))
          itmp=a%i
          jtmp=a%j
          vtmp=a%val
          deallocate(a%i,a%j,a%val)
          allocate(a%i(nnz+1),a%j(nnz+1),a%val(nnz+1))
          a%i(1:nnz)=itmp
          a%j(1:nnz)=jtmp
          a%val(1:nnz)=vtmp
          a%i(nnz+1)=i
          a%j(nnz+1)=j
          a%val(nnz+1)=val
          deallocate(itmp,jtmp,vtmp)
        end if
      end if
    end subroutine setcoeff





    ! fonction d'ajout d'une valeur a un coefficient
    !     si la matrice est vide, elle contient cet unique coefficient
    subroutine addtocoeff(a, i, j, val)
      implicit none
      integer, intent(in) :: i, j
      type(matsparse), intent(inout) :: a
      real(kind=8), intent(in) :: val
      integer, dimension(:), pointer :: itmp
      integer, dimension(:), pointer :: jtmp
      real(kind=8), dimension(:), pointer :: vtmp
      integer :: nnz, ind
      logical :: isEmpty

      if (i>a%n) then
        stop "ERROR in setcoeff: row index out of bounds" 
      end if
      if (j>a%m) then
        stop "ERROR in setcoeff: column index out of bounds"
      end if

      isEmpty=.false.
      if (a%isAllocated.eqv..false.) then
        call msallocate(a,1)
        isEmpty=.true.
      end if

      if (size(a%i) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays i and j in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%i) /= size(a%val)) then
        stop "ERROR in setcoeff: arrays i and val in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%val) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays j and val in matsparse must have the same size (did you use msallocate ?)"
      end if

      if (isEmpty.eqv..true.) then
        a%i(1)=i
        a%j(1)=j
        a%val(1)=val
      else
        ind=findpos(a, i, j)
        if (ind /= 0) then ! coeff already exists
          a%val(ind)=a%val(ind)+val
        else ! coeff does not exists, we add it to the storage
          nnz=size(a%i)
          allocate(itmp(nnz), jtmp(nnz), vtmp(nnz))
          itmp=a%i
          jtmp=a%j
          vtmp=a%val
          deallocate(a%i,a%j,a%val)
          allocate(a%i(nnz+1),a%j(nnz+1),a%val(nnz+1))
          a%i(1:nnz)=itmp
          a%j(1:nnz)=jtmp
          a%val(1:nnz)=vtmp
          a%i(nnz+1)=i
          a%j(nnz+1)=j
          a%val(nnz+1)=val
          deallocate(itmp,jtmp,vtmp)
        end if
      end if
    end subroutine addtocoeff




    ! fonction de multiplication d'un coefficient par une valeur
    subroutine multcoeff(a, i, j, val)
      implicit none
      integer, intent(in) :: i, j
      type(matsparse), intent(inout) :: a
      real(kind=8), intent(in) :: val
      integer, dimension(:), pointer :: itmp
      integer, dimension(:), pointer :: jtmp
      real(kind=8), dimension(:), pointer :: vtmp
      integer :: nnz, ind
      logical :: isEmpty

      if (i>a%n) then
        stop "ERROR in setcoeff: row index out of bounds" 
      end if
      if (j>a%m) then
        stop "ERROR in setcoeff: column index out of bounds"
      end if

      isEmpty=.false.
      if (a%isAllocated.eqv..false.) then
        call msallocate(a,1)
        isEmpty=.true.
      end if

      if (size(a%i) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays i and j in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%i) /= size(a%val)) then
        stop "ERROR in setcoeff: arrays i and val in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%val) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays j and val in matsparse must have the same size (did you use msallocate ?)"
      end if

      if (isEmpty.eqv..true.) then
        a%i(1)=i
        a%j(1)=j
        a%val(1)=val
      else
        ind=findpos(a, i, j)
        if (ind /= 0) then ! coeff already exists
          a%val(ind)=a%val(ind)*val
        else ! coeff does not exists, we add it to the storage
          nnz=size(a%i)
          allocate(itmp(nnz), jtmp(nnz), vtmp(nnz))
          itmp=a%i
          jtmp=a%j
          vtmp=a%val
          deallocate(a%i,a%j,a%val)
          allocate(a%i(nnz+1),a%j(nnz+1),a%val(nnz+1))
          a%i(1:nnz)=itmp
          a%j(1:nnz)=jtmp
          a%val(1:nnz)=vtmp
          a%i(nnz+1)=i
          a%j(nnz+1)=j
          a%val(nnz+1)=val
          deallocate(itmp,jtmp,vtmp)
        end if
      end if
    end subroutine multcoeff




    ! fonction de suppresion d'un coefficient
    subroutine delcoeff(a, i, j)
      implicit none
      integer, intent(in) :: i, j
      type(matsparse), intent(inout) :: a
      integer, dimension(:), pointer :: itmp
      integer, dimension(:), pointer :: jtmp
      real(kind=8), dimension(:), pointer :: vtmp
      integer :: nnz, ind

      if (i>a%n) then
        stop "ERROR in setcoeff: row index out of bounds"
      end if
      if (j>a%m) then
        stop "ERROR in setcoeff: column index out of bounds"
      end if

      if (a%isAllocated.eqv..false.) then
        return
      end if

      if (size(a%i) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays i and j in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%i) /= size(a%val)) then
        stop "ERROR in setcoeff: arrays i and val in matsparse must have the same size (did you use msallocate ?)"
      end if
      if (size(a%val) /= size(a%j)) then
        stop "ERROR in setcoeff: arrays j and val in matsparse must have the same size (did you use msallocate ?)"
      end if

      ind=findpos(a, i, j)
      if (ind /= 0) then ! coeff exists

        nnz=size(a%i)
        allocate(itmp(nnz-1), jtmp(nnz-1), vtmp(nnz-1))
        itmp(1:ind-1)=a%i(1:ind-1)
        itmp(ind:nnz-1)=a%i(ind+1:nnz)
        jtmp(1:ind-1)=a%j(1:ind-1)
        jtmp(ind:nnz-1)=a%j(ind+1:nnz)
        vtmp(1:ind-1)=a%val(1:ind-1)
        vtmp(ind:nnz-1)=a%val(ind+1:nnz)
        deallocate(a%i,a%j,a%val)
        allocate(a%i(nnz-1),a%j(nnz-1),a%val(nnz-1))
        a%i=itmp
        a%j=jtmp
        a%val=vtmp
        deallocate(itmp,jtmp,vtmp)
      end if
    end subroutine delcoeff



    ! fonction de recuperation d'un coefficient
    function coeff(a, i, j) result(val)
      type(matsparse), intent(in) :: a
      integer, intent(in) :: i, j
      real(kind=8) :: val
      integer :: ind

      ind=findpos(a, i, j)
      if (ind /= 0) then
        val=a%val(ind)
      else
        val=0.d0
      end if
    end function coeff




    ! fonction de definition de la taille de la matrice
    subroutine setsize(a,n,m)
      type(matsparse), intent(inout) :: a
      integer, intent(in) :: n, m
      a%n=n
      a%m=m
    end subroutine setsize




    ! fonction d'affichage d'une matrice
    subroutine affiche(a)
      implicit none
      type(matsparse), intent(in) :: a
      integer :: i
      print *, "matrice de taille", a%n, a%m
      print *, "nnz", size(a%i)
      do i=1,size(a%i)
        print *,a%i(i),a%j(i),a%val(i)
      end do
    end subroutine affiche
    ! fonction de copie d'une matrice
    !     la matrice resultat doit etre vide
    subroutine spcopy(a,b)
      implicit none
      type(matsparse), intent(out) :: a
      type(matsparse), intent(in) :: b
      integer :: nnz

      if (b%isAllocated.eqv..true.) then
        nnz=size(b%i)
        call sparse(a,b%n,b%m,nnz)
        a%i=b%i
        a%j=b%j
        a%val=b%val
      else
        call sparse(a,b%n,b%m)
      end if
      a%isAllocated=b%isAllocated
    end subroutine spcopy




    ! fonction d'affectation de valeur a tous les coefficients d'une matrice
    subroutine spaffect(a,b)
      implicit none
      type(matsparse), intent(out) :: a
      real(kind=8), intent(in) :: b
      a%val=b
    end subroutine spaffect
    ! recherche la position d'un coefficient
    function findpos(a, i, j, mo) result(ind)
      implicit none
      type(matsparse), intent(in) :: a
      integer, intent(in) :: i, j
      integer, intent(in), optional :: mo
      integer :: ind, nnz, k, m

      nnz=size(a%i)
      m=nnz
      if (present(mo)) then
        m=mo
      end if
      ind=0
      do k=1,m
        if ((a%i(k) == i).and.(a%j(k) == j)) then
          ind=k
          exit
        end if
      end do
    end function findpos




    ! recherche des coefficients correspondant a un critere
    function findind(a, mask) result(ind)
      implicit none
      logical, dimension(:), intent(in) :: mask
      type(matsparse), intent(in) :: a
      integer, dimension(:), pointer :: ind
      integer :: nnz, k, m

      nnz=size(a%i)
      if (nnz /= size(mask)) then
        stop "Erreur, dimension incompatible du mask"
      end if

      m=count(mask)
      allocate(ind(m))
      ind=pack((/(k,k=1,nnz)/),mask)
    end function findind




    ! extrait les coefficients ayant une certaine valeur
    ! permet par exemple de ne recuperer que les coefficients reellement non nuls
    function extract(a,mask) result(b)
      implicit none
      type(matsparse), intent(in) :: a
      logical, dimension(:), intent(in) :: mask
      type(matsparse) :: b
      integer :: nnza, nnzb

      nnza=size(a%i)
      if (nnza /= size(mask)) then
        stop "Erreur, dimension incompatible du mask"
      end if
      nnzb=count(mask)
      call sparse(b,a%n,a%m,nnzb)
      b%i=pack(a%i,mask)
      b%j=pack(a%j,mask)
      b%val=pack(a%val,mask)
    end function extract




    ! tri d'une matrice
    subroutine sort(a, sens)
      implicit none
      type(matsparse), intent(inout) :: a
      integer, optional, intent(in) :: sens
      real(kind=8) :: tmp
      integer :: itmp, iglob, jglob, k, choix
      logical :: ok
      choix=1
      if (present(sens).and.((sens==0).or.(sens==1))) then
        choix=sens
      end if
      do ! Tri
        ok = .true.
        do k=2,size(a%i)
          ! calcul des indices globaux
          if (choix==1) then
            iglob=(a%i(k-1)-1)*a%n+a%j(k-1)
            jglob=(a%i(k)-1)*a%n+a%j(k)
          else
            iglob=(a%j(k-1)-1)*a%m+a%i(k-1)
            jglob=(a%j(k)-1)*a%m+a%i(k)
          end if
          if (iglob > jglob) then
            ok = .false.
            ! echange des indices de ligne
            itmp=a%i(k-1)
            a%i(k-1)=a%i(k)
            a%i(k)=itmp
            ! echange des indices de colonne
            itmp=a%j(k-1)
            a%j(k-1)=a%j(k)
            a%j(k)=itmp
            ! echange des valeurs
            tmp=a%val(k-1)
            a%val(k-1)=a%val(k)
            a%val(k)=tmp
          end if
        end do
        if (ok) exit
      end do
    end subroutine sort




    ! transposition d'une matrice
    function sptranspose(a) result(b)
      type(matsparse), intent(in) :: a
      type(matsparse) :: b

      call sparse(b,a%m,a%n,size(a%i))
      b%i=a%j
      b%j=a%i
      b%val=a%val
      call sort(b)
    end function sptranspose




    ! addition de deux matrices (standard)
    function spadd(a,b) result(c)
      implicit none
      type(matsparse), intent(in) :: a, b
      type(matsparse) :: c, tmp
      integer :: k, l, nnz, nnzmax
      integer, dimension(size(b%i)) :: done, ind
      integer, dimension(:), allocatable :: indfinal
      logical :: found

      if ((a%n /= b%n).or.(a%m /= b%m)) then
        stop "dimensions de matrices incompatibles"
      end if

      done=1
      ind=(/(k,k=1,size(b%i))/)
      nnzmax=size(a%i)+size(b%i)

      call sparse(tmp,a%n,a%m,nnzmax)
      k=1
      nnz=1
      do k=1,size(a%i)
        found=.false.
        l=1
        do while(l<=size(b%i) .and. (found .eqv. .false.))
          if((a%i(k) == b%i(l)) .and. (a%j(k) == b%j(l))) then
            if((a%val(k)+b%val(l)) /= 0.d0) then
              tmp%i(nnz)=a%i(k)
              tmp%j(nnz)=a%j(k)
              tmp%val(nnz)=a%val(k)+b%val(l)
              nnz=nnz+1
            end if
            found=.true.
            done(l)=0
          end if
          l=l+1
        end do
        if (found .eqv. .false.) then
          tmp%i(nnz)=a%i(k)
          tmp%j(nnz)=a%j(k)
          tmp%val(nnz)=a%val(k)
          nnz=nnz+1
        end if
      end do

      k=sum(done)
      if(k/=0) then
        allocate(indfinal(k))
        indfinal=pack(ind,done==1)
        do l=1,k
          tmp%i(nnz)=b%i(indfinal(l))
          tmp%j(nnz)=b%j(indfinal(l))
          tmp%val(nnz)=b%val(indfinal(l))
          nnz=nnz+1
        end do
        deallocate(indfinal)
      end if
      call sparse(c,a%n,a%m,nnz-1)
      c%i=tmp%i(1:nnz-1)
      c%j=tmp%j(1:nnz-1)
      c%val=tmp%val(1:nnz-1)
      deallocate(tmp%i,tmp%j,tmp%val)
    end function spadd





    ! addition de deux matrices (par tri et fusion)
    function spadd2(a,b) result (c)
      implicit none
      type(matsparse), intent(in) :: a, b
      type(matsparse) :: c, tmp
      integer :: k, l, nnza,nnzb, nnzmax, nnz, iglob, jglob
      logical, dimension(:), allocatable ::mask

      if ((a%n /= b%n).or.(a%m /= b%m)) then
        stop "dimensions de matrices incompatibles"
      end if

      nnza=size(a%i)
      nnzb=size(b%i)
      nnzmax=nnza+nnzb
      ! print *, nnza, nnzb, nnzmax
      allocate(tmp%i(nnzmax),tmp%j(nnzmax),tmp%val(nnzmax),mask(nnzmax))
      tmp%n=a%n
      tmp%m=a%m
      tmp%i(1:nnza)=a%i
      tmp%i(nnza+1:nnzmax)=b%i
      tmp%j(1:nnza)=a%j
      tmp%j(nnza+1:nnzmax)=b%j
      tmp%val(1:nnza)=a%val
      tmp%val(nnza+1:nnzmax)=b%val
      call sort(tmp)

      do k=2,nnzmax
        iglob=(tmp%i(k-1)-1)*tmp%n + tmp%j(k-1)
        jglob=(tmp%i(k)-1)*tmp%n + tmp%j(k)
        if (iglob == jglob) then
          tmp%val(k-1)=tmp%val(k-1)+tmp%val(k)
          tmp%val(k)=0.d0
        end if
      end do

      mask=tmp%val/=0.d0
      nnz=count(mask)
      allocate(c%i(nnz),c%j(nnz),c%val(nnz))
      c%n=a%n
      c%m=a%m
      c%i=pack(tmp%i,mask)
      c%j=pack(tmp%j,mask)
      c%val=pack(tmp%val,mask)
      deallocate(mask,tmp%i,tmp%j,tmp%val)
    end function spadd2




    ! addition de deux matrices pre-triees
    function spaddsort(a,b) result(c)
      implicit none
      type(matsparse), intent(in) :: a, b
      type(matsparse) :: c, tmp
      integer :: k, l, nnz, nnzmax,nnza, nnzb, i

      if ((a%n /= b%n).or.(a%m /= b%m)) then
        stop "dimensions de matrices incompatibles"
      end if
      nnza=size(a%i)
      nnzb=size(b%i)
      nnzmax=nnza+nnzb

      allocate(tmp%i(nnzmax),tmp%j(nnzmax),tmp%val(nnzmax))
      k=1
      l=1
      nnz=1
      do while ((k <= nnza) .and. (l <= nnzb))
        if (a%i(k) < b%i(l)) then
          tmp%i(nnz)=a%i(k)
          tmp%j(nnz)=a%j(k)
          tmp%val(nnz)=a%val(k)
          k=k+1
          nnz=nnz+1
        else if (a%i(k) == b%i(l)) then
          if (a%j(k) < b%j(l)) then
            tmp%i(nnz)=a%i(k)
            tmp%j(nnz)=a%j(k)
            tmp%val(nnz)=a%val(k)
            k=k+1
            nnz=nnz+1
          else if (a%j(k) == b%j(l)) then
            if ((a%val(k) + b%val(l)) /= 0.d0) then
              tmp%i(nnz)=a%i(k)
              tmp%j(nnz)=a%j(k)
              tmp%val(nnz)=a%val(k) + b%val(l)
              nnz=nnz+1
            end if
            k=k+1
            l=l+1
          else
          tmp%i(nnz)=b%i(l)
          tmp%j(nnz)=b%j(l)
          tmp%val(nnz)=b%val(l)
          l=l+1
          nnz=nnz+1
          end if
        else
          tmp%i(nnz)=b%i(l)
          tmp%j(nnz)=b%j(l)
          tmp%val(nnz)=b%val(l)
          l=l+1
          nnz=nnz+1
        end if
      end do
      ! traitement final derniers elements si positionne a un endroit different
      if(l <= nnzb) then
        do i=l,nnzb
          tmp%i(nnz)=b%i(i)
          tmp%j(nnz)=b%j(i)
          tmp%val(nnz)=b%val(i)
          nnz=nnz+1
        end do
      end if
      if (k <= nnza) then
        do i=k,nnza
          tmp%i(nnz)=a%i(i)
          tmp%j(nnz)=a%j(i)
          tmp%val(nnz)=a%val(i)
          nnz=nnz+1
        end do
      end if

      allocate(c%i(nnz-1),c%j(nnz-1),c%val(nnz-1))
      c%n=a%n
      c%m=a%m
      c%i=tmp%i(1:nnz-1)
      c%j=tmp%j(1:nnz-1)
      c%val=tmp%val(1:nnz-1)
      deallocate(tmp%i,tmp%j,tmp%val)
    end function spaddsort




    ! soustraction de deux matrices
    function spminus(a,b) result (c)
      implicit none
      type(matsparse), intent(in) :: a, b
      type(matsparse) :: c, tmp
      integer :: k, l, nnza,nnzb, nnzmax, nnz, iglob, jglob
      logical, dimension(:), allocatable ::mask

      if ((a%n /= b%n).or.(a%m /= b%m)) then
        stop "dimensions de matrices incompatibles"
      end if

      nnza=size(a%i)
      nnzb=size(b%i)
      nnzmax=nnza+nnzb

      call sparse(tmp,a%n,a%m,nnzmax)
      tmp%i(1:nnza)=a%i
      tmp%i(nnza+1:nnzmax)=b%i
      tmp%j(1:nnza)=a%j
      tmp%j(nnza+1:nnzmax)=b%j
      tmp%val(1:nnza)=a%val
      tmp%val(nnza+1:nnzmax)=-b%val
      call sort(tmp)

      do k=2,nnzmax
        iglob=(tmp%i(k-1)-1)*tmp%n + tmp%j(k-1)
        jglob=(tmp%i(k)-1)*tmp%n + tmp%j(k)
        if (iglob == jglob) then
          tmp%val(k-1)=tmp%val(k-1) + tmp%val(k)
          tmp%val(k)=0.d0
        end if
      end do

      mask=tmp%val/=0.d0
      nnz=count(mask)
      call sparse(c,a%n,a%m,nnz)
      c%i=pack(tmp%i,mask)
      c%j=pack(tmp%j,mask)
      c%val=pack(tmp%val,mask)
      deallocate(mask,tmp%i,tmp%j,tmp%val)
    end function spminus




    ! effectue le produit d'une matrice par un scalaire
    function spscalmat(a,b) result(c)
      implicit none
      type(matsparse), intent(in) :: a
      type(matsparse) :: c
      real(kind=8), intent(in) :: b
      integer :: nnz
      if(b==0.d0) then
        stop "multiplication d'une sparse par 0"
      end if
      nnz=size(a%i)

      call sparse(c,a%n,a%m,nnz)
      c%i=a%i
      c%j=a%j
      c%val=b*a%val
    end function spscalmat
    ! effectue le produit d'un scalaire par une matrice
    function spmatscal(b,a) result(c)
      implicit none
      type(matsparse), intent(in) :: a
      type(matsparse) :: c
      real(kind=8), intent(in) :: b
      integer :: nnz

      if(b==0.d0) then
        stop "multiplication d'une sparse par 0"
      end if
      nnz=size(a%i)

      call sparse(c,a%n,a%m,nnz)
      c%i=a%i
      c%j=a%j
      c%val=b*a%val
    end function spmatscal




    ! effectue le produit d'une matrice par un vecteur
    function spmatvec(a,b) result(u)
      implicit none
      type(matsparse), intent(in) :: a
      real(kind=8), dimension(:), intent(in) :: b
      real(kind=8), dimension(a%n) :: u
      integer :: i,nnz

      if (a%m /= size(b)) then
        stop "dimensions matrice / vecteur incompatibles"
      end if
      u=0.d0
      nnz=size(a%i)
      do i=1,nnz
        u(a%i(i))=u(a%i(i))+a%val(i)*b(a%j(i))
      end do
    end function spmatvec




    ! effectue le produit d'une matrice par une matrice
    function spmatmat(a,b) result(u)
      implicit none
      type(matsparse), intent(in) :: a, b
      type(matsparse) :: u,tmp
      integer :: i,j,k,nnz,nnzmax,nnza,nnzb

      if (a%m /= b%n) then
        stop "dimensions de matrices incompatibles"
      end if
      nnza=size(a%i)
      nnzb=size(b%i)
      nnzmax=nnza*nnzb
      call sparse(tmp,a%n,a%m,nnzmax)
      tmp%val=0.d0
      nnz=1
      do i=1,nnza
        do j=1,nnzb
          if (a%j(i) == b%i(j)) then
            k=find(tmp,a%i(i),b%j(j),nnz-1)
            if(k/=0) then
              tmp%i(k)=a%i(i)
              tmp%j(k)=b%j(j)
              tmp%val(k)=tmp%val(k)+a%val(i)*b%val(j)
            else
              tmp%i(nnz)=a%i(i)
              tmp%j(nnz)=b%j(j)
              tmp%val(nnz)=tmp%val(nnz)+a%val(i)*b%val(j)
              nnz=nnz+1
            end if
          end if
        end do
      end do
      u=extract(tmp,tmp%val/=0.d0)
      deallocate(tmp%i,tmp%j,tmp%val)
    end function spmatmat




    ! lufact - calcule factorisation LU: A=L*U
    subroutine lufact(A, L, U, klo, kuo)
      implicit none
      type(matsparse), intent(in) :: A
      type(matsparse), intent(out) :: L, U
      type(matsparse) :: Ltmp, Utmp
      integer, optional, intent(in) :: klo, kuo
      integer :: i, j, k, nu, nl, n, nnzl, nnzu, kl, ku, ind1, ind2

      n=A%n

      kl=A%n-1
      if(present(klo).and.klo<n) then
        kl=klo
      end if

      ku=A%n-1
      if(present(kuo).and.kuo<n) then
        ku=kuo
      end if
      nnzl=(kl+1)*(2*n-kl)/2
      nnzu=(ku+1)*(2*n-ku)/2
      call sparse(Ltmp, A%n, A%m, nnzl)
      call sparse(Utmp, A%n, A%m, nnzu)
      nnzl=1
      nnzu=1

      do k=1,n
        ! traitement du coefficient diagonal
        Ltmp%i(nnzl)=k
        Ltmp%j(nnzl)=k
        Ltmp%val(nnzl)=1.d0
        nnzl=nnzl+1
        ind1=find(A,k,k)
        if (ind1/=0) then
          Utmp%val(nnzu)=A%val(ind1)
        else
          Utmp%val(nnzu)=0.d0
        end if
        do j=1,k-1
          ind1=find(Ltmp,k,j,nnzl)
          ind2=find(Utmp,j,k,nnzu)
          if (ind1/=0.and.ind2/=0) then
            Utmp%val(nnzu)=Utmp%val(nnzu)-Ltmp%val(ind1)*Utmp%val(ind2)
          end if
        end do
        if(Utmp%val(nnzu)/=0.d0) then
          Utmp%i(nnzu)=k
          Utmp%j(nnzu)=k
          nnzu=nnzu+1
        end if

        ! traitement des coefficients sous-diagonaux de L
        nl=min(k+kl,n)
        do i=k+1,nl
          ind1=find(A,i,k)
          if(ind1/=0) then
            Ltmp%val(nnzl)=A%val(ind1)
          end if
          do j=1,k-1
            ind1=find(Ltmp,i,j,nnzl)
            ind2=find(Utmp,j,k,nnzu)

            if (ind1/=0.and.ind2/=0) then
              Ltmp%val(nnzl)=Ltmp%val(nnzl)-Ltmp%val(ind1)*Utmp%val(ind2)
            end if
          end do
          if(Ltmp%val(nnzl)/=0.d0) then
            ind1=find(Utmp,k,k,nnzu)
            Ltmp%i(nnzl)=i
            Ltmp%j(nnzl)=k
            Ltmp%val(nnzl)=Ltmp%val(nnzl)/Utmp%val(ind1)
            nnzl=nnzl+1
          end if
        end do

        ! traitement des coefficients sur-diagonaux de U
        nu=min(k+ku,n)
        do i=k+1,nu
          ind1=find(A,k,i)
          if(ind1/=0) then
            Utmp%val(nnzu)=A%val(ind1)
          else
            Utmp%val(nnzu)=0.d0
          end if
          do j=1,k-1
            ind1=find(Ltmp,k,j,nnzl)
            ind2=find(Utmp,j,i,nnzu)
            if (ind1/=0.and.ind2/=0) then
              Utmp%val(nnzu)=Utmp%val(nnzu)-Ltmp%val(ind1)*Utmp%val(ind2)
            end if
          end do
          if(Utmp%val(nnzu)/=0.d0) then
            Utmp%i(nnzu)=k
            Utmp%j(nnzu)=i
            nnzu=nnzu+1
          end if
        end do
      end do
      L=Ltmp
      U=Utmp
      ! numeriquement un coeff nul est un coeff plus petit,
      ! en valeur absolue, que le zero machine
      ! L=extract(Ltmp,dabs(Ltmp%val)>epsilon(1.d0))
      ! U=extract(Utmp,dabs(Utmp%val)>epsilon(1.d0))
      deallocate(Ltmp%i,Ltmp%j,Ltmp%val)
      deallocate(Utmp%i,Utmp%j,Utmp%val)
    end subroutine lufact




    ! lusolve - resout le systeme L*U x = f par L*y=f puis U*x=y
    subroutine lusolve(L,U,f,x,klo,kuo)
      implicit none
      type(matsparse), intent(in) :: L, U
      real(kind=8), dimension(:), intent(in) :: f
      real(kind=8), dimension(size(f)), intent(out) :: x
      real(kind=8), dimension(size(f)) :: y
      integer, optional, intent(in) :: klo, kuo
      integer :: i, j, n, kl, ku, nl, nu, ind

      n=size(f)

      kl=size(f)-1
      if(present(klo).and.klo<n) then
        kl=klo
      end if

      ku=size(f)-1
      if(present(kuo).and.kuo<n) then
        ku=kuo
      end if

      x=0.d0
      y=0.d0

      ! resolution de L*y=f
      do i=1,n
        y(i)=f(i)
        nl=max(1,i-kl)
        do j=nl,i-1
          ind=find(L,i,j)
          if(ind/=0) then
            y(i)=y(i)-L%val(ind)*y(j)
          end if
        end do
      end do

      ! resolution de U x = y
      do i=n,1,-1
        x(i)=y(i)
        nu=min(n,i+ku)
        do j=i+1,nu
          ind=find(U,i,j)
          if (ind/=0) then
            x(i)=x(i)-U%val(ind)*x(j)
          end if
        end do
        ind=find(U,i,i)
        x(i)=x(i)/U%val(ind)
      end do
    end subroutine lusolve
    ! definition de la matrice identite sparse
    function speye(n) result(a)
      implicit none
      integer, intent(in) :: n
      type(matsparse) :: a
      integer :: i

      allocate(a%i(n),a%j(n),a%val(n))
      a%n=n
      a%m=n
      a%i=(/(i,i=1,n)/)
      a%j=(/(i,i=1,n)/)
      a%val=1.d0
    end function speye
    ! routine de test


    ! fonction de descente pour une matrice triangulaire inférieure
    function downSolve(L,y) result(x)

      implicit none

      type(matsparse), intent(in)                        :: L
      real(kind=8), dimension(:), intent(in)             :: y
      real(kind=8), dimension(size(y))                   :: x
      integer :: i, kl, nl, n, ind, j

      ! taille pour la résolution
      n = size(y)

      ! resolution de L*x=y
      do i=1,n
        x(i)=y(i)
        do j=1,i-1
          ind=find(L,i,j)
          if(ind/=0) then
            x(i)=x(i)-L%val(ind)*x(j)
          end if
        end do
        ind = find(L, i, i)
        if( L%val(ind) /= 0) then
           x(i) = x(i)/L%val(ind)
        else
           write(*,*) "ERROR : la matrice n'est pas inversible"
        end if
     end do

    end function downSolve


    subroutine test_amsta01sparse()
      type(matsparse) :: mat, mat2, mat3, mat4, mat5, mat6, mat7, mat8, L, U, res
      integer, parameter :: n=10
      real(kind=8) :: eps
      integer, dimension(:), allocatable :: ind, indth
      real(kind=8), dimension(:), allocatable :: x,xth, b
      integer :: i
      ! recherche du zero machine
      eps=epsilon(1.d0)

      print *,
      print *, "-------------- TEST MODULE AMSTA01SPARSE -----------"
      print *,

      print *, "------------------------------------------"
      print *, "Test setcoeff"

      call sparse(mat,2,2)
      call setcoeff(mat,1,1,2.d0)
      call sparse(res,2,2,1)
      res%i(1)=1
      res%j(1)=1
      res%val(1)=2.d0
      if (dsqrt(sum((mat%i-res%i)*(mat%i-res%i)+(mat%j-res%j)*(mat%j-res%j)+(mat%val-res%val)*(mat%val-res%val))) < eps) then
        print *, "     test setcoeff - first coeff OK"
      else
        print *, "     erreur setcoeff - first coeff"
      end if

      call setcoeff(mat,1,2,-1.d0)
      deallocate(res%i,res%j,res%val)
      call sparse(res,2,2,2)
      res%i=(/1,1/)
      res%j=(/1,2/)
      res%val=(/2.d0,-1.d0/)

      if (dsqrt(sum((mat%i-res%i)*(mat%i-res%i)+(mat%j-res%j)*(mat%j-res%j)+(mat%val-res%val)*(mat%val-res%val))) < eps) then
        print *, "     test setcoeff - insertion OK"
      else
        print *, "     error setcoeff - insertion"
      end if

      call setcoeff(mat,1,2,1.d0)
      res%i=(/1,1/)
      res%j=(/1,2/)
      res%val=(/2.d0,1.d0/)
      if (dsqrt(sum((mat%i-res%i)*(mat%i-res%i)+(mat%j-res%j)*(mat%j-res%j)+(mat%val-res%val)*(mat%val-res%val))) < eps) then
        print *, "     test setcoeff - update OK"
      else
        print *, "     error setcoeff - update"
      end if

      print *, "------------------------------------------"
      print *, "Test addition"
      call msreallocate(mat,5)
      call setsize(mat,3,3)
      mat%i=(/ 2, 1, 2, 2, 3/)
      mat%j=(/ 1, 1, 3, 2, 1/)
      mat%val=(/ 2.d0, 1.d0, -2.d0, -1.d0, 1.d0/)
      call sparse(mat2,3,3,6)
      mat2%i=(/ 2, 1, 3, 1, 3, 3/)
      mat2%j=(/ 1, 1, 1, 2, 3, 2/)
      mat2%val=(/ -2.d0, 2.d0, 2.d0, 3.d0, -1.d0, 3.d0/)

      mat3=spadd(mat,mat2)
      call msreallocate(res,7)
      res%i=(/ 1, 2, 2, 3, 1, 3, 3 /)
      res%j=(/ 1, 3, 2, 1, 2, 3, 2 /)
      res%val=(/ 3.d0, -2.d0, -1.d0, 3.d0, 3.d0, -1.d0, 3.d0/)

      if(dsqrt(sum((mat3%i-res%i)*(mat3%i-res%i)+(mat3%j-res%j)*(mat3%j-res%j)+(mat3%val-res%val)*(mat3%val-res%val))) < eps) then
        print *, "     test addition - unsorted OK"
      else
        print *, "     error addition - unsorted"
      end if

      mat4=spadd2(mat,mat2)
      res%i=(/ 1, 1, 2, 2, 3, 3, 3 /)
      res%j=(/ 1, 2, 2, 3, 1, 2, 3 /)
      res%val=(/ 3.d0, 3.d0, -1.d0, -2.d0, 3.d0, 3.d0, -1.d0/)

      if(dsqrt(sum((mat4%i-res%i)*(mat4%i-res%i)+(mat4%j-res%j)*(mat4%j-res%j)+(mat4%val-res%val)*(mat4%val-res%val))) < eps) then
        print *, "     test addition - sort and merge OK"
      else
        print *, "     error addition - sort and merge"
      end if

      call sort(mat)
      call sort(mat2)
      mat5=spaddsort(mat2,mat)
      res%i=(/ 1, 1, 2, 2, 3, 3, 3 /)
      res%j=(/ 1, 2, 2, 3, 1, 2, 3 /)
      res%val=(/ 3.d0, 3.d0, -1.d0, -2.d0, 3.d0, 3.d0, -1.d0/)

      if(dsqrt(sum((mat5%i-res%i)*(mat5%i-res%i)+(mat5%j-res%j)*(mat5%j-res%j)+(mat5%val-res%val)*(mat5%val-res%val))) < eps) then
        print *, "     test addition - sorted OK"
      else
        print *, "     error addition - sorted"
      end if

      print *, "------------------------------------------"
      print *, "Test findind"
      ind=find(mat5,mat5%val==3.d0)
      allocate(indth(4))
      indth=(/ 1, 2, 5, 6 /)
      if(size(ind)==size(indth) .and. dsqrt(dble(sum((indth-ind)*(indth-ind)))) < eps) then
        print *, "     test findind OK"
      else
        print *, "     erreur findind"
      end if

      print *, "------------------------------------------"
      print *, "Test findpos"
      i=find(mat5,2,3)
      if (i==4) then
        print "('      test findpos OK : i=2, j=3, imax=7, k=',i1,'=4')",i
      else
        print "('      erreur findpos : i=2, j=3, imax=7, k=',i1,'/=4')",i
      end if
      i=find(mat5,2,3,3)
      if (i==0) then
        print "('      test findpos OK : i=2, j=3, imax=3, k=',i1,'=0')",i
      else
        print "('      erreur findpos : i=2, j=3, imax=3, k=',i1,'/=0')",i
      end if
      i=find(mat5,1,3)
      if (i==0) then
        print "('      test findpos OK : i=1, j=3, imax=7, k=',i1,'=0')",i
      else
        print "('      erreur findpos : i=1, j=3, imax=7, k=',i1,'/=0')",i
      end if

      print *, "------------------------------------------"
      print *,"Test lufact"
      call msreallocate(mat,4)
      call msreallocate(mat2,5)

      mat%i=(/ 1, 2, 2, 3 /)
      mat%j=(/ 1, 1, 2, 3 /)
      mat%val=(/ 1.d0, 2.d0, 1.d0, 1.d0 /)
      mat2%i=(/ 1, 1, 2, 2, 3 /)
      mat2%j=(/ 1, 3, 2, 3, 3 /)
      mat2%val=(/ -1.d0, 1.d0, -2.d0, 1.d0, -1.d0 /)

      mat6=mat*mat2
      call lufact(mat6,L,U)
      mat7=(-1.d0)*(L*U)
      mat8=mat7+mat6

      if( dsqrt(sum((mat8%val)*(mat8%val))) > eps) then
        print *, "     error lufact"
      else
        print *, "     test lufact OK"
      end if

      print *, "------------------------------------------"
      print *, "Test lusolve"
      allocate(x(3),xth(3),b(3))
      xth=(/ 1.d0, 2.d0, 3.d0 /)
      b=(/ 2.d0, 3.d0, -3.d0 /)
      call lusolve(L, U, b, x)

      if(dsqrt(sum((x-xth)*(x-xth))) > eps) then
        print *, "     error lusolve"
      else
        print *, "     test lusolve OK"
      end if

      deallocate(mat%i,mat%j,mat%val)
      deallocate(mat2%i,mat2%j,mat2%val)
      deallocate(mat3%i,mat3%j,mat3%val)
      deallocate(mat4%i,mat4%j,mat4%val)
      deallocate(mat5%i,mat5%j,mat5%val)
      deallocate(mat6%i,mat6%j,mat6%val)
      deallocate(mat7%i,mat7%j,mat7%val)
      deallocate(mat8%i,mat8%j,mat8%val)
      deallocate(L%i,L%j,L%val)
      deallocate(U%i,U%j,U%val)
      deallocate(ind)
    end subroutine test_amsta01sparse
end module amsta01sparse
