!-----------------------------------------------------------------------------
! module de résolution du problème de laplacien par éléments avec Jacobi et Gauss-Sield
! Auteur : Emeriau PE
!-----------------------------------------------------------------------------

module amsta01solveur

  use amsta01maillage
  use amsta01sparse
  use amsta01probleme

  implicit none

contains


  ! *************************************************************************
  !                            JACOBI
  ! on veut resoudre AX = B avec A = M - N où A diag et N à diag vide
  ! *************************************************************************

  subroutine solveJacobi(pb, eps, conv, nbSsDomaines, myRank, ierr)

    !include 'mpif.h'

    !implicit none

    ! variables d'entrée et de sortie du problème
    type(probleme), intent(inout)          :: pb           ! problème que l'on veut résoudre
    real, intent(in)                       :: eps          ! critère de précision sur la convergence
    integer, intent(in)                    :: myRank, ierr
    logical, intent(out)                   :: conv         ! variable logique pour tester la convergence
    integer, intent(in)                    :: nbSsDomaines ! Nombre de ss-domaines

    ! variables locales
    type(matsparse)                       :: N, M_inv                    ! avec K=M-N
    real(kind=8), dimension(:), pointer   :: uk, rk                      ! itéré kieme de la solution et du résidu
    real(kind=8), dimension(:), pointer   :: uk_prime, uk_sec            ! contient les noeuds à envoyer pour les communications
    real(kind=8)                          :: carre_norm, norm, norm_init ! norme du résidu
    integer                               :: n_size, i, k, j             ! taille vecteur solution, variables boucles

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
       !else
       !   call setcoeff(M_inv, i, i, 0d0)
       !end if
       end if
    end do

    ! Initialisation du vecteur solution
    uk = 0.d0

    ! Vecteur residu initial
    rk = spmatvec(pb%p_Kelim, uk) - pb%felim

    ! Produit scalaire pour un domaine donne
    carre_norm = dot_product(rk,rk)

    ! On fait la somme et on redistribue a tout le monde
    call MPI_ALLREDUCE(carre_norm, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! norme intiale
    norm_init = dsqrt(norm)

    
    
    ! on alloue le vecteur uk_prime qui contient les noeuds à envoyer
    allocate(uk_prime(size(pb%mesh%int2glob)))
    ! On alloue le vecteur uk_sec
    if (myRank /= 0) allocate(uk_sec(size(pb%mesh%intFront2glob)))
    if (myRank == 0) allocate(uk_sec(size(pb%mesh%intFront2glob_proc0(1,:))))


    ! Boucle :
    do k = 1, 5000 ! Pour ne pas avoir de boucle infinie (on pourrait optimiser ce critère sachant que k <= n_size ?)

       ! calcul de l'itéré kieme de la solution
       uk = M_inv * (N * uk) + M_inv * pb%felim    ! spmatvec -> *

       ! *********************************************************
       ! *********************************************************
       ! COMMUNICATIONS PARALLELLES :

       ! *** 1)  Interface vers ss-domaines :


       ! on remplit le vecteur uk_prime des noeuds à envoyer grâce à int2glob sur le proc 0 (interface)
       if (myRank == 0) uk_prime(:) = uk(pb%mesh%int2glob(:))
       ! on envoit uk_prime de 0 vers les autres proc
       call MPI_Bcast(uk_prime, size(uk_prime), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       ! on réaffecte chaque noeud de l'interface pour uk sur chaque processeur
       uk(pb%mesh%int2glob(:)) = uk_prime(:)


       ! *** 2) Ss-domaines vers interface :


       ! Envoi et réception des noeuds voisins de la frontière grâce aux vecteurs intFront2glob
       ! et le tableau récapitulatif intFront2glob_proc0
       if (myRank == 0) then
          do i = 1, nbSsDomaines
             uk_sec = 0
             call MPI_RECV(uk_sec, size(uk_sec), MPI_DOUBLE_PRECISION, i, 100, MPI_COMM_WORLD, status, ierr)

             do j = 1, size(pb%mesh%intFront2glob_proc0(i,:))
                if(pb%mesh%intFront2glob_proc0(i,j) /= 0) then
                   uk(pb%mesh%intFront2glob_proc0(i,j)) = uk_sec(j)
                end if
             end do

          end do
       else
          uk_sec(:) = uk(pb%mesh%intFront2glob(:))
          call MPI_SEND(uk_sec, size(uk_sec), MPI_DOUBLE_PRECISION, 0, 100, MPI_COMM_WORLD, ierr)
       end if


       ! *** 3) Calcul du résidu
       rk = pb%felim - pb%p_Kelim * uk         ! spmatvec -> *

       ! Si le sommet n'appartient pas au sous domaine alors on met rk a 0
       ! On pourrait avoir des problemes pour des sommets aux frontieres
       carre_norm = dot_product(rk,rk)
       call MPI_ALLREDUCE(carre_norm, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       norm = dsqrt(norm)

       ! sortie de la boucle si on a atteint la convergence
       if(norm < eps*norm_init) then
          conv = .TRUE.
          if (myRank == 0) then
             write(*,*) 'INFO    : Residu reel : ', norm
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
    call MPI_REDUCE(uk, pb%u, n_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! désallocation des matrices crées
    deallocate(uk, rk)

  end subroutine solveJacobi










  ! **********************************************************************************
  !                            GAUSS-SEIDEL
  ! on veut resoudre AX = B avec A = M - N où M triang inf et N triang sup à diag vide
  ! **********************************************************************************
  
  subroutine solveGaussSeidel(pb, eps, conv, myRank, nbSsDomaine, ierr)

    implicit none

    type(probleme), intent(inout) :: pb
    real, intent(in)              :: eps
    logical, intent(out)          :: conv
    integer, intent(in)           :: myRank, ierr  ! Variables MPI
    integer, intent(in)           :: nbSsDomaine

    ! Variable locale MPI
    integer, dimension(MPI_STATUS_SIZE)         :: status

    ! Variables locales
    type(matsparse)                      :: M, N, Add, AdDp
    real(kind=8), dimension(:), pointer  :: rk, uk
    real(kind=8), dimension(:), pointer  :: uk_prime, uk_sec, uk_tri
    real(kind=8)                         :: norm, norm_init, sum    
    integer                              :: n_size,i,k,j, errcode
    logical :: supp



    ! conv est mise a false par default
    conv = .FALSE.
    supp = .TRUE.

    ! Osn recupere la taille du probleme avec elimination
    n_size = size(pb%felim)

    ! On alloue les valeurs des vecteurs solution et residu
    allocate(uk(n_size), rk(n_size))


    call spcopy(Add,  pb%p_Kelim)
    call spcopy(AdDp, pb%p_Kelim) 


    boucle_chCol0 : do j=1,n_size

       do i=1,n_size
          if (coeff(Add,j,i) /= 0) supp = .FALSE.
       end do

       if (supp .eqv. .TRUE.) then 
          do i=1,n_size
             call delcoeff(Add,i,j)
          end do
       else
          do i=1,n_size
             call delcoeff(AdDp,i,j)
          end do
       end if

       supp = .TRUE.

    end do boucle_chCol0

    ! Definition des matrices M et N. Attention K = M - N !
    N = spmatscal(-1.d0, extract(Add, Add%i < Add%j))
    M = extract(Add, Add%i >= Add%j)



    ! On enleve les donnees inutiles dans felim
    do j=1,n_size
       if(pb%mesh%refPartNodes(j) /= myRank) pb%felim(j) = 0
    end do


    ! Initialisation du vecteur solution
    uk = 0.d0

    ! Vecteur residu initial
    rk = spmatvec(pb%p_Kelim, uk) - pb%felim

    ! Produit scalaire pour un domaine donne
    sum = dot_product(rk,rk)

    ! On fait la somme et on redistribue a tout le monde
    call MPI_ALLREDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! norme intiale
    norm_init = dsqrt(norm)


    ! On alloue le vecteur uk_prime qui contient les noeuds a envoyer
    allocate(uk_prime(size(pb%mesh%int2glob)))

    ! On alloue le vecteur uk_sec qui contient les noeuds a envoyer
    if(myRank /= 0) allocate(uk_sec(size(pb%mesh%intFront2glob)))
    if(myRank == 0) allocate(uk_sec(size(pb%mesh%intFront2glob_proc0(1,:))))



    ! On preferera faire une boucle do pour ne pas avoir de fuite. On sort avec un exit.
    do  k = 1,5000

       if (myRank /= 0) then

          ! Iteration de uk
          uk = spmatvec(N,uk) + pb%felim - spmatvec(AdDp,uk) 
          uk = downSolve(M,uk,.TRUE.)


          uk_sec(:) = uk(pb%mesh%intFront2glob(:))
          call MPI_SEND(uk_sec, size(pb%mesh%intFront2glob(:)), &
               MPI_DOUBLE_PRECISION, 0, 100, MPI_COMM_WORLD, ierr)

       else if (myRank == 0) then 

          uk = spmatvec(N,uk) + pb%felim

          do i=1,nbSsDomaine

             uk_sec = 0

             call MPI_RECV(uk_sec, size(uk_sec), &
                  MPI_DOUBLE_PRECISION, i, 100, MPI_COMM_WORLD, status, ierr)

             do j = 1,size(pb%mesh%intFront2glob_proc0(i,:))
                if(pb%mesh%intFront2glob_proc0(i,j) /= 0) uk(pb%mesh%intFront2glob_proc0(i,j)) = uk_sec(j)
             end do

          end do

          ! Alorithme de descente modifié dans notre cas
          uk = downsolve(M, uk - spmatvec(AdDp,uk),.TRUE., uk, AdDp)

       end if


       ! On remplit le vecteur uk_prime des noeuds a envoyer grace à int2glob sur le proc 0 (interface)
       if(myRank == 0) uk_prime(:) = uk(pb%mesh%int2glob(:))
       ! on envoit uk_prime de 0 vers les autres proc
       call MPI_BCAST(uk_prime, size(uk_prime), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       ! on reaffecte chaque noeud de l'interface pour uk sur chaque processeur
       if(myRank /= 0) uk(pb%mesh%int2glob(:)) = uk_prime(:)



       ! Calcul de la norme de du residu pour observer la convergence
       ! On fait ce calcul toutes les 10 iterations pour aller plus vite. Utile ?
       if (mod(k,10) == 0) then

          ! Calcul de residu et de la norme
          rk = spmatvec(pb%p_Kelim, uk) - pb%felim

          ! Produit scalaire pour un domaine donne
          sum = dot_product(rk,rk)

          ! On fait la somme et on redistribue a tout le monde
          call MPI_ALLREDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

          ! Calcul de la norme
          norm = dsqrt(norm)

          ! Si jamais on a atteint le critère de convergence on sort de la boucle
          if (norm < eps*norm_init) then
             conv = .TRUE.
             if(myRank == 0) write(*,*) 'INFO    : Residu reel : ', norm
             if(myRank == 0) write(*,*) 'INFO    : Precision attendue pour la convergence : ', eps
             if(myRank == 0) write(*,*) 'INFO    : Convergence apres ', k, ' iterations de la methode de Gauss Seidel'
             exit
          end if

       end if
    end do


    ! On recompose la solution avec les donnees de chaque processeur
    ! Si le sommet n'appartient pas au sous domaine alors on met uk a 0
    ! On pourrait avoir des problemes pour des sommets aux frontieres
    do j=1,n_size
       if(pb%mesh%refPartNodes(j) /= myRank) uk(j) = 0
    end do

    ! On recupere tout sur le processeur 0
    call MPI_REDUCE(uk, pb%u, n_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! On desallocate les matrice creees
    deallocate(uk,rk,uk_prime,uk_sec)

  end subroutine solveGaussSeidel



end module amsta01solveur
