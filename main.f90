program main

  use amsta01maillage
  use amsta01sparse
  use amsta01probleme
  use amsta01solveur

  use mpi

  implicit none


  ! variables relatives à la définition du problème
  type(maillage)                       :: mail
  type(probleme)                       :: pb
  type(matsparse)                      :: Kt, Mt
  real(kind=8)                         :: erreur
  real(kind=8), dimension(:), pointer  :: residu
  logical                              :: conv   ! Indique s'il y a eu convergence pour les methodes iterative
  integer                              :: nbSsDomains
  character(len=20)                    :: filename

  !variables relatives à l'utilisation de mpi
  integer                                :: nbTask, myRank, ierr, request, errcode
  integer, dimension(MPI_STATUS_SIZE)    :: status



  ! initialisation MPI
  call MPI_Init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  write(*,*)
  write(*,*) '  **** TA01 Equation de la chaleur ****  '

  open(unit=11, file="python_res.txt", form='formatted')

  read(11,*) filename
  read(11,*) nbSsDomains

  write(*,*)
  write(*,*) '_________________________________________'
  write(*,*) 'Nombre de sous-domaines du maillage lu grâce au script Python :'
  write(*,*) nbSsDomains

  ! erreur si le nombre de sous-domaines est différent de celui du nombre de processeurs
  if(nbTask /= nbSsDomains + 1) then
     if(myRank == 0) then
        write(*,*) "___________________________________________________________________________________"
        write(*,*) "ERROR : Le nombre de sous-domaines est différent du nombre de processeurs demandés"
        write(*,*) "Le programme va s'arrêter"
        write(*,*) "___________________________________________________________________________________"
     end if
     call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
  end if



  write(*,*)
  write(*,*) '_________________________________________'
  write(*,*) 'Proprietes du maillage :'

  ! lecture du maillage
  mail = loadFromMshFile("./testpart.msh", nbSsDomains)

  ! construction des donnees sur les triangles
  call getTriangles(mail,  myRank, nbSsDomains)


  ! Affichage des données des noeuds et des elements
  if(myRank == 0) call affichePartNoeud(mail, "infoNoeuds.log")
  if(myRank == 0) call affichePartElem(mail, "infoElems.log")
  call affichePartTri(mail, "infoTris.log", myRank)

  ! creation du probleme
  call loadFromMesh(pb,mail)

  ! assemblage des matrices elements finis
  call assemblage(pb)

  ! pseudo-elimination des conditions essentielles
  call pelim(pb,mail%refNodes(1))

  write(*,*) '_________________________________________'
  write(*,*) 'Erreur theorique attendu :'

  ! calcul du residu theorique
  allocate(residu(mail%nbNodes))
  residu=pb%felim-pb%p_Kelim*pb%uexa
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Erreur theorique=", erreur


  write(*,*) '_________________________________________'
  write(*,*) 'Resolution du systeme lineaire : '

  ! Resolution par jacobi
  call solveJacobi(pb, 0.000001, conv)

  ! Resolution par Gauss Seidel
  !call solveGaussSeidel(pb, 0.000001, conv)

  ! Si on n'a pas converge on utilise une methode directe
  if (conv .eqv. .FALSE.) then
     ! resolution du systeme lineaire
     call solveLU(pb)
     write(*,*) 'WARNING : Il n y a pas eu convergence de la methode iterative'
     write(*,*) 'INFO    : Le systeme a ete resolu a l aide d une methode directe LU'
  end if


  write(*,*) '_________________________________________'
  write(*,*) 'Calcul du residu reel et de l erreur :'

  ! calcul du residu
  residu=pb%felim-pb%p_Kelim*pb%u
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Residu=", erreur

  ! calcul de l'erreur L2
  erreur=dsqrt(dot_product(pb%uexa-pb%u,pb%uexa-pb%u))
  print *, "||u-uexa||_2=", erreur

  write(*,*) '_________________________________________'
  write(*,*)
  write(*,*) '      **** Fin du programmme ****'
  write(*,*)

  ! sauvegarde de la solution et de la solution theorique
  call saveToVtu(pb%mesh,pb%u,pb%uexa)


call MPI_Finalize(ierr)

end program
