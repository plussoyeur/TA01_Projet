!-----------------------------------------------------------------------------
! module de résolution du problème de laplacien par éléments avec Jacobi et Gauss-Sield
! Auteur : Emeriau PE
!-----------------------------------------------------------------------------

module amsta01solveur

  use amsta01maillage
  use amsta01sparse
  use amsta01probleme

  contains


  ! calcul de la solution du probleme par Jacobi
  ! on veut resoudre AX = B avec A = M - N où A diag et N tradiag à diag vide
  subroutine solveJacobi(pb, eps, conv, nbSsDomaines, myRank, ierr)

    !include 'mpif.h'
    
    implicit none

    ! variables d'entrée et de sortie du problème
    type(probleme), intent(inout)          :: pb           ! problème que l'on veut résoudre
    real, intent(in)                       :: eps          ! critère de précision sur la convergence
    integer, intent(in)                    :: myRank, ierr
    logical, intent(out)                   :: conv         ! variable logique pour tester la convergence
    integer, intent(in)                    :: nbSsDomaines ! Nombre de ss-domaines
    
    ! variables locales
    type(matsparse)                       :: N, M_inv                  ! avec K=M-N
    real(kind=8), dimension(:), pointer   :: uk, rk                    ! itéré kieme de la solution et du résidu
    real(kind=8), dimension(:), pointer   :: uk_prime, uk_sec, uk_tri  ! contient les noeuds à envoyer pour les communications
    real(kind=8)                          :: carre_norm, norm          ! norme du résidu
    integer                               :: n_size, i, k, j           ! taille vecteur solution, variables boucles

    ! variable locale MPI
    integer, dimension(MPI_STATUS_SIZE)   :: status

    ! on définit la convergence à faux au départ
    conv = .FALSE.

    ! on prend la taille du problème à résoudre (avec élimination)
    n_size = size(pb%felim)

    ! on alloue les vecteurs itérés au rang k de la solution et du résidu
    allocate(uk(n_size), rk(n_size))

    ! Definition des matrices M_inv et N
    N = (-1.d0) * extract(pb%p_Kelim, pb%p_Kelim%i /= pb%p_Kelim%j)    ! Attention on a K = M - N (spcopy -> =; spscalmat -> *)
    call sparse(M_inv, n_size, n_size)
    ! Remplissage et suppresion de M_inv
    ! Lorsque le coefficient de la matrice M_inv n'est pas déf (car le noeud n'appartient pas
    ! au processeur), on le met à 0 pour que le produit matriciel réalisé ensuite se passe
    ! sans problème (on met certaines valeurs en trop sur les proc des ss-domaines mais multiplier
    ! par 0 en faisant cela.
    do i = 1, n_size
       if(coeff(pb%p_Kelim,i,i) /= 0) then
          call setcoeff(M_inv, i, i,  (1.d0)/(coeff(pb%p_Kelim, i, i)))
       else
          call setcoeff(M_inv, i, i, 0d0)
       end if
    end do

    allocate(uk_tri(n_size))
    ! on alloue le vecteur uk_prime qui contient les noeuds à envoyer
    allocate(uk_prime(size(pb%mesh%int2glob)))
    ! On alloue le vecteur uk_sec 
    if (myRank /= 0) allocate(uk_sec(size(pb%mesh%intFront2glob)))
    if (myRank == 0) allocate(uk_sec(size(pb%mesh%intFront2glob_proc0(1,:))))

    
    ! Initialisation du vecteur solution et boucle
    uk = 1.d0
    do k = 1, 5000 ! Pour ne pas avoir de boucle infinie (on pourrait optimiser ce critère sachant que k <= n_size ?)

       ! calcul de l'itéré kieme de la solution
       uk = M_inv * (N * uk) + M_inv * pb%felim    ! spmatvec -> *

       uk_tri = uk
       ! *********************************************************
       ! *********************************************************
       ! COMMUNICATIONS PARALLELLES :

       ! *** 1)  Interface vers ss-domaines :


       ! on remplit le vecteur uk_prime des noeuds à envoyer grâce à int2glob sur le proc 0 (interface)
       if (myRank == 0) uk_prime(:) = uk(pb%mesh%int2glob(:))
       ! on envoit uk_prime de 0 vers les autres proc
       call MPI_Bcast(uk_prime, size(uk_prime), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
       ! on réaffecte chaque noeud de l'interface pour uk sur chaque processeur
       uk(pb%mesh%int2glob(:)) = uk_prime(:)


       ! *** 2) Ss-domaines vers interface : 


       ! Envoi et réception des noeuds voisins de la frontière grâce aux vecteurs intFront2glob
       ! et le tableau récapitulatif intFront2glob_proc0
       if (myRank == 0) then
          do i = 1, nbSsDomaines
             uk_sec = 0
             call MPI_RECV(uk_sec, size(uk_sec), MPI_DOUBLE, i, 100, MPI_COMM_WORLD, status, ierr)

             do j = 1, size(pb%mesh%intFront2glob_proc0(i,:))
                if(pb%mesh%intFront2glob_proc0(i,j) /= 0) then
                   uk(pb%mesh%intFront2glob_proc0(i,j)) = uk_sec(j)
                end if
             end do
             
          end do
       else
          uk_sec(:) = uk(pb%mesh%intFront2glob(:))
          call MPI_SEND(uk_sec, size(uk_sec), MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, ierr)
       end if
       
       
       ! *********************************************************
       ! *********************************************************

       if (myRank == 0 .AND. k == 5) then
          do j=1,n_size
             if (uk(j) /= uk_tri(j)) then
                write(*,*) ' Changement pour : ', j
             end if
          end do
       end if
         
       
       ! *** 3) Calcul du résidu
       rk = pb%felim - pb%p_Kelim * uk         ! spmatvec -> *

       ! Si le sommet n'appartient pas au sous domaine alors on met rk a 0
       ! On pourrait avoir des problemes pour des sommets aux frontieres
       do j=1,n_size
          if(pb%mesh%refPartNodes(j) /= myRank) rk(j) = 0
       end do

       carre_norm = dot_product(rk,rk)
       call MPI_ALLREDUCE(carre_norm, norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
       norm = dsqrt(norm)
       
       ! sortie de la boucle si on a atteint la convergence
       if(norm < eps) then
          conv = .TRUE.
          if (myRank == 0) then 
             write(*,*)
             write(*,*) 'INFO    : Precision attendue pour la convergence : ', eps
             write(*,*) 'INFO    : Convergence apres ', k, ' iterations de la methode de Jacobi'
          end if
          exit
       end if

    end do

    ! On recompose la solution avec les données de chaque processeur
    ! Si le sommet n'appartient pas au sous domaine alors on met le noeud correspondant
    ! de uk à 0
    do j = 1, n_size
       if(pb%mesh%refPartNodes(j) /= myRank) uk(j) = 0
    end do
    
    ! On récupère tout sur le processeur 0
    call MPI_REDUCE(uk, pb%u, n_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! désallocation des matrices crées
    deallocate(uk, rk)

  end subroutine solveJacobi




  ! calcul de la solution du probleme par Gauss Seidel
  ! on veut resoudre AX = B avec A = M - N où A triang inf et N triang strict sup
  subroutine solveGaussSeidel(pb, eps, conv, myRank)

    ! variables d'entrée du problème et de sortie
    type(probleme), intent(inout)     :: pb     ! probleme que l'on veut resoudre
    real, intent(in)                  :: eps    ! critere de convergence
    integer, intent(in)               :: myRank
    logical, intent(out)              :: conv   ! variable logique pour tester la convergence
   

    ! variables locales
    type(matsparse)                     :: N, M_inv          ! avec K=M-N
    integer                             :: n_size, k, i      ! taille du vecteur solution, entiers pour les boucles
    real(kind=8), dimension(:), pointer :: uk, rk            ! vecteur solution à l'itération k et résidu à l'ordre k
    real(kind=8)                        :: norm              ! norme du résidu à l'ordre k


    ! initialisation à faux de la variable logique de convergence
    conv = .FALSE.

    ! définition de la taille du problème (avec élimination) à résoudre et allocation des vecteurs uk et rk
    n_size = size(pb%felim)
    allocate(uk(n_size), rk(n_size))

    ! définition des matrices M_inv et N
    call sparse(M_inv, n_size, n_size)
    call sparse(N, n_size, n_size)

    N = (-1.d0) * extract(pb%p_Kelim, pb%p_Kelim%i < pb%p_Kelim%j)
    M_inv = extract(pb%p_Kelim, pb%p_Kelim%i >= pb%p_Kelim%j)

    ! algo de descente pour le calcul de uk
    ! Initialisation du vecteur solution et boucle
    uk = 1.d0
    do k = 1, 1000 ! Pour ne pas avoir de boucle infinie (on pourrait optimiser ce critère sachant que k <= n_size ?)

       ! calcul de l'itéré kieme de la solution
       uk = (N * uk) + pb%felim    ! spmatvec -> *
       uk = downSolve(M_inv, uk)
       ! calcul du résidu
       rk = pb%felim - pb%p_Kelim * uk         ! spmatvec -> *
       norm = dsqrt(dot_product(rk,rk))

       ! sortie de la boucle si on a atteint la convergence
       if(norm < eps) then
          conv = .TRUE.
          if (myRank == 0) then
             write(*,*)
             write(*,*) 'INFO    : Precision attendue pour la convergence : ', eps
             write(*,*) 'INFO    : Convergence apres ', k, ' iterations de la methode de Gauss-Seidel'
          end if
          exit
       end if

    end do

    
   ! la solution est la valeur du dernier itéré
    pb%u = uk

  end subroutine solveGaussSeidel









end module amsta01solveur
