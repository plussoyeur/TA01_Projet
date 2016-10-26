! ----------------------------------------------------------------
! module de lecture de maillages gmsh
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------
module amsta01maillage

  implicit none

  ! classe maillage
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! refElems contient le nombre de sous-domaines auxquels touchent l'élément (ex : 2 si appartient à 1 et touche 3)
  !
  ! elemsPartRef contient les infos realtives aux sous-domaines auxquels touchent l'élément (ex : 1 -3 si appartient à 1 et touche 3)
  !
  ! RefPartNodes associe le numéro du noeud avec son sous-domaine et permet donc de tester si un numéro de noeud est associé à plusieurs sous-domaines -> noeud appartenant à l'interface
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  type maillage
    integer :: nbNodes, nbElems, nbTri
    real(kind=8), dimension(:,:), pointer :: coords
    integer, dimension(:,:), pointer :: typeElems, elemsVertices, triVertices, elemsPartRef, RefPartTri
    integer, dimension(:), pointer :: refNodes, refElems, refTri, elemsNbPart, triNbPart, RefPartNodes
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
         write(*,*) "boucle while"
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
    subroutine getTriangles(mail, nbSsDomains)

      implicit none

      type(maillage), intent(inout) :: mail
      integer, intent(in)           :: nbSsDomains
      integer                       :: i, j, nbTri
      nbTri=0

      do i=1, mail%nbElems
         if (mail%typeElems(i,1) == 2) then
            nbTri=nbTri+1
         end if
      end do
      mail%nbTri=nbTri
      Print*, "NbTri=", NbTri

      allocate(mail%refTri(nbTri), mail%triVertices(nbTri,3), mail%RefPartTri(nbTri,nbSsDomains), mail%triNbPart(nbTri))
      j=1
      do i=1, mail%nbElems
         if (mail%typeElems(i,1) == 2) then
            mail%refTri(j)=mail%refElems(i)
            mail%triVertices(j,1:3)=mail%elemsVertices(i,1:3)
            mail%triNbPart(j) = mail%elemsNbPart(i)
            mail%RefPartTri(j,1:mail%triNbPart(j)) = mail%elemsPartRef(i,1:mail%elemsNbPart(i))
            j=j+1
         end if
      end do
    end subroutine getTriangles


  ! construction d'une subroutine de stockage des éléments des tableaux
  subroutine affichePart(mail, filename)

    implicit none

    type(maillage), intent(in)     :: mail
    character(len=*), intent(in)   :: filename
    integer                        :: j

    open(unit=17, file=filename, form='formatted')

    ! stockage dans un fichier des informations associées aux noeuds
    do j=1,mail%nbNodes
       write(17,*) 'Noeud numero : ', j, ' | RefNodes : ', mail%refNodes(j), ' | RefPartNodes : ', mail%refPartNodes(j)
    end do

    ! affichage des informations associées aux triangles
    !do j=1,mail%nbTri
    !     write(*,*) 'Noeud numero : ', j, ' | RefNodes : ', mail%refTri(j), ' | RefPartNodes : ', mail%RefPartTri(j,:)
    !end do


    close(17)

  end subroutine affichePart


end module amsta01maillage
