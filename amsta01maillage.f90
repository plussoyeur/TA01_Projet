! ----------------------------------------------------------------
! module de lecture de maillages gmsh
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------
module amsta01maillage

  implicit none

  type maillage
    integer :: nbNodes, nbElems, nbTri
    real(kind=8), dimension(:,:), pointer :: coords
    integer, dimension(:,:), pointer :: typeElems, elemsVertices, triVertices, elemsPartRef, RefPartTri
    integer, dimension(:), pointer :: refNodes, refElems, refTri, elemsNbPart, triNbPart, RefPartNodes, tri2elem
  end type


  contains


    ! retourne le nombre de noeuds d'un type de maille gmsh donne
    function type2nbNodes(t) result(tt)
      implicit none
      integer, intent(in) :: t
      integer :: tt
      integer, dimension(31) :: nbNodesOfTypes
      nbNodesOfTypes = (/2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13,9,10,12,15,15,21,4,5,6,20,35,56 /)
      tt=nbNodesOfTypes(t)
    end function



    ! construit un maillage par lecture d'un fichier gmsh
    function loadFromMshFile(filename, nbSsDomains) result(res)
      implicit none

      ! déclaration des variables
      character(len=*), intent(in)   :: filename
      integer, intent(in)            :: nbSsDomains
      type(maillage)                 :: res
      integer                        :: ios, ibuf, ibufD, ibufT, i, j, nbtags
      character(len=100)             :: sbuf, sbuf2, sbuf3
      real(kind=8)                   :: rbuf
      integer, dimension(:), pointer :: elemData

      ! ouverture du fichier et test de son existence
      open(unit=10, file=filename, form='formatted', status='old', iostat=ios)
      if (ios /= 0) then
         stop "ERROR File not opened"
      end if

      do while (ios == 0)
         ! lecture de la première ligne
         read(10,*, iostat=ios) sbuf


         if (sbuf == "$MeshFormat") then
            ! read next line containing the file format version, and 2 integers about numeric precision
            read(10,*, iostat=ios) sbuf, sbuf2, sbuf3
            ! check if block ends correctly
            read(10,'(A14)', iostat=ios) sbuf

            if (sbuf /= "$EndMeshFormat") then
               stop "ERROR : Msh $MeshFormat block must end with $EndMeshFormat"
            end if


         else if (sbuf =="$Nodes") then
            ! read number of nodes
            read(10,*, iostat=ios) res%nbNodes
            print*, "nbNodes=",res%nbNodes
            allocate(res%coords(res%nbNodes,3), res%refNodes(res%nbNodes))
            res%refNodes=0

            ! read coordinates
            do i=1,res%nbNodes
               read(10,*, iostat=ios) ibuf, res%coords(i,1), res%coords(i,2), res%coords(i,3)
            end do

            ! check if block ends correctly
            read(10,'(A9)', iostat=ios) sbuf

            if (sbuf /= "$EndNodes") then
               stop "ERROR : Msh $Nodes block must end with $EndNodes"
            end if


         else if (sbuf =="$Elements") then
            read(10,*, iostat=ios) res%nbElems
            print*, "nbElems=",res%nbElems
            allocate(res%typeElems(res%nbElems,2), res%refElems(res%nbElems), res%elemsVertices(res%nbElems,3))
            allocate(res%elemsPartRef(res%nbElems,nbSsDomains), res%refPartNodes(res%nbNodes), res%elemsNbPart(res%nbElems))

            ! initialisation de RefPartNodes à 0 (! hors de la boucle)
            res%refPartNodes = 0

            ! read elements data
            do i=1,res%nbElems
               ! We get each line in a character string
               read(10,'(A)', iostat=ios) sbuf
               ! convert character string to array of integers
               allocate(elemData(100))
               read(sbuf,*,iostat=ios) elemData
               ! get type of element and number of vertices related to it
               res%typeElems(i,1)=elemData(2)
               res%typeElems(i,2)=type2nbNodes(res%typeElems(i,1))

               ! association éléments / sous-domaines auxquels il appartient
               res%elemsNbPart(i) = elemdata(6)
               do j = 1, res%elemsNbPart(i)
                  res%elemsPartRef(i,j) = elemData(6+j)
               end do


               ! management of tags (domain ref of element and eventually partition data)
               nbtags=elemData(3)
               res%refElems(i)=elemData(4)

               ! read vertices of element
               do j=1,res%typeElems(i,2)
                  res%elemsVertices(i,j)=elemData(nbtags+3+j)
               end do


               ! set reference of vertices
               do j=1,res%typeElems(i,2)

                  if (res%refNodes(res%elemsVertices(i,j)) == 0) then
                     res%refNodes(res%elemsVertices(i,j)) = res%refElems(i)
                  end if

               end do

               ! affectation des noeuds à un domaine et teste de leur appartenance à plusieurs domaines
               ! 1 -> noeud du bord; 2 -> noeud de l'intérieur; -7 -> noeud d'une interface; -3 -> noeud interface inter bord
               do j = 1,res%typeElems(i,2)

                  ! on teste si le noeud a déjà été lu et on l'attribue à un sous-domaine
                  if(res%RefPartNodes(res%elemsVertices(i,j)) == 0) then
                     res%RefPartNodes(res%elemsVertices(i,j)) = elemData(7)

                  ! on teste maintenant si le domaine lu est différent de celui enregistré
                  else if(res%RefPartNodes(res%elemsVertices(i,j)) /= elemData(7)) then
                     ! S'il s'agit d'un noeud sur le bord on le met à -3
                     if(res%refNodes(res%elemsVertices(i,j)) == 1) res%refNodes(res%elemsVertices(i,j)) = -3
                     ! S'il s'agit d'un noeud du bord inter interface on le met à -7
                     if (res%refNodes(res%elemsVertices(i,j)) >= 0) res%refNodes(res%elemsVertices(i,j)) = -7

                  end if
               end do

               deallocate(elemData)

            end do

            ! check if block ends correctly
            read(10,'(A12)', iostat=ios) sbuf

            if (sbuf /= "$EndElements") then
               stop "ERROR : Msh $Elements block must end with $EndElements"
            end if

         else

         end if

      end do

      do j = 1, res%nbNodes

         if(res%refNodes(j) < 0) res%refPartNodes(j) = 0

      end do

      close(10)

    end function




















    ! construit la liste des triangles du maillage
    subroutine getTriangles(mail, myRank, nbSsDomaine)

      implicit none

      type(maillage), intent(inout)   :: mail
      integer, optional, intent(in)   :: nbSsDomaine
      integer, intent(in)             :: myRank
      integer                         :: i, j, k, nbTri_tot

      nbTri_tot   = count(mail%typeElems(:,1) == 2)
      mail%nbTri  = nbTri_tot
      if(myRank == 0) Print*, "NbTri (total) =", mail%nbTri

      allocate(mail%refTri(mail%nbTri), mail%refPartTri(mail%nbTri,nbSsDomaine))
      allocate(mail%triNbPart(mail%nbTri), mail%tri2elem(mail%nbTri))

      ! Initialisation du tableau de partition des triangles a 0
      mail%refPartTri = 0
      j = 1

      ! Identification des triangles parmis les elements
      boucle_identification_triangle : do i=1,mail%nbElems

         condition_est_un_triangle : if (mail%typeElems(i,1) == 2) then

            mail%triNbPart(j) = mail%elemsNbpart(i)
            mail%refTri(j)=mail%refElems(i)
            mail%refPartTri(j,1:mail%triNbPart(j)) = mail%elemsPartRef(i,1:mail%elemsNbPart(i))
            mail%tri2elem(j) = i
            j=j+1

         end if condition_est_un_triangle

      end do boucle_identification_triangle

      !! Note ---------
      ! Le tableau tri2elem est un tableau qui permet lorsque l'on connait l'identifiant d'un
      ! triangle de remonter a l'identifiant de l'element associe. Dans la boucle precedente
      ! le triangle j correspond a l'element numerote i. Ceci nous fait stocker un tableau de
      ! plus mais ce n'est qu'un vecteur et il simplifie les choses considerablement ensuite.
      !! --------------


      ! On change la variable nbTri pour quelle corresponde a ce qu'elle vaudra
      ! non plus au global mais sur chaque processeur consideres et on alloue
      ! le tableau des identifiants des sommets des triangles du processeur
      ! en fonction de cette variable

      mail%nbTri = count(mail%refPartTri(:,1) == myRank)
      allocate(mail%triVertices(mail%nbTri,3))

      j = 1

      ! Récupération des identifiants des somments
      boucle_triVertices_proc_myRank : do i=1,nbTri_tot
         if (mail%refPartTri(i,1) == myRank) then
            mail%triVertices(j,1:3) = mail%elemsVertices(mail%tri2elem(i),1:3)
            j=j+1
         end if
      end do boucle_triVertices_proc_myRank



      condition_rank_0 : if (myRank == 0) then

         ! Initialisation du nombre de triangles pour le processeur 0 a 0
         mail%nbTri = 0

         ! Récupération du nombre de triangles touchant l'interface
         do i=1,nbTri_tot
            ! On compte le nombre de triangle ayant au moins un noeud touchant le bord
            mail%nbTri = mail%nbTri + &
                 min(count(mail%refPartNodes(mail%elemsVertices(mail%tri2elem(i),:))==0),1)
         end do

         ! Allocation de la matrice triVertices pour le processeur 0
         allocate(mail%triVertices(mail%nbTri,3))

         k = 1

         ! Attribution des numéros des sommets en fonction des numerosnuméros des éléments
         do i = 1,nbTri_tot
            do j=1,3
               ! On regarde si l'un des noeuds associe a un refPartNodes de 0 (i.e. il est sur l'interface)
               if (mail%refPartNodes(mail%elemsVertices(mail%tri2elem(i),j))==0) then
                  mail%triVertices(k,1:3) = mail%elemsVertices(mail%tri2elem(i),1:3)
                  k = k+1
                  ! Une fois que le triangle a été attribué à cause de l'un de ces sommets on sort
                  ! ceci pour ne pas compter deux fois un meme triangle
                  exit
               end if
            end do
         end do

      end if condition_rank_0


    end subroutine getTriangles
















    !! -------------------------------------------------------------- !!
    !! Subroutines d'affichage
    !! -------------------------------------------------------------- !!



    ! Affiche les références des noeuds
    subroutine affichePartNoeud(mail, filename)

      implicit none

      type(maillage), intent(in) :: mail
      character(len=*), intent(in) :: filename
      integer    :: j


      open(unit=13, file=filename, form="formatted")

      write(13,*) '!!!! -------- Info Noeuds -------- !!!!'
      write(13,*)
      write(13,*) 'Il y a ', mail%nbNodes, ' noeuds'
      write(13,*)

      do j=1,mail%nbNodes
         write(13,*) 'Noeud numero : ', j, ' | RefNodes : ', &
              mail%refNodes(j), ' |  RefPartNodes : ', mail%refPartNodes(j)
      end do

      close(13)


    end subroutine affichePartNoeud






    ! Affiche les références des éléments
    subroutine affichePartElem(mail, filename)

      implicit none

      type(maillage), intent(in) :: mail
      character(len=*), intent(in) :: filename
      integer    :: j


      open(unit=12, file=filename, form="formatted")

      write(12,*) '!!!! -------- Info elems -------- !!!!'
      write(12,*)
      write(12,*) 'Il y a ', mail%nbElems, ' elems'
      write(12,*)

      do j=1,mail%nbElems
         write(12,*) 'Element numero : ', j, ' | RefElems : ', &
              mail%refElems(j), ' |  RefPartElems : ', mail%elemsPartRef(j,:), &
              '| Type d element : ', mail%typeElems(j,1)
      end do

      close(12)


    end subroutine affichePartElem






    ! Affiche les références des triangles
    subroutine affichePartTri(mail, filename, myRank)

      implicit none

      type(maillage), intent(in)     :: mail
      character(len=*), intent(in)   :: filename
      integer, intent(in)            :: myRank
      character(len=100)             :: filename_bis, filename_tris
      character(len=2)               :: str_rank
      integer, dimension(2)          :: nb_tri
      integer    :: j, k

      write(str_rank, '(I1.1)') myRank

      k = index(filename, '.log')
      filename_bis  = trim(filename(1:k-1)//'_proc_'//trim(str_rank)//'.log')
      filename_tris = trim(filename(1:k-1)//'_total_.log')


      open(unit=14, file=filename_bis, form="formatted")

      write(14,*) '!!!! -------- Info tris  -------- !!!!'
      write(14,*)
      write(14,*) 'Il y a ', mail%nbTri, ' triangles pour le proc ', myRank
      write(14,*)

      do j=1,mail%nbTri
         write(14,*) 'Triangle ''numero'' : ', j, ' | triVertices : ', &
              mail%triVertices(j,:)
      end do

      close(14)


      if (myRank == 0) then

         open(unit=15, file=filename_tris, form="formatted")

         nb_tri = shape(mail%refPartTri)

         write(15,*) '!!!! -------- Info tris  -------- !!!!'
         write(15,*)
         write(15,*) 'Il y a ', nb_tri(1), ' triangles au total'
         write(15,*)

         do j=1,nb_tri(1)
            write(15,*) 'Triangle numero : ', mail%tri2elem(j), ' | refPartTri : ', &
                 mail%refPartTri(j,:)
         end do

         close(15)

      end if


    end subroutine affichePartTri




end module amsta01maillage
