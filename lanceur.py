# -*-coding:utf-8 -*

import os
import sys
import numpy


if len(sys.argv) > 1:
    filename = sys.argv[1]
    
else:

    filename = raw_input("Entrer le nom du fichier mesh \n")

    # Si jamais filename n'est pas specifie on lui donne testpart.msh
    # Valeur par defaut pour gagner du temps sur les tests
    if filename  == "":
        filename = "testpart.msh"



        
# Test de l'existence du fichier de maillage
try:
    with open(filename): pass
except IOError:
    print "*************************************************************"
    print "      ERREUR : le fichier spécifié n'existe pas              "
    print "*************************************************************"
    sys.exit(1)



file = open(filename, "r")

while 1:
    f = file.readline()
    if list(f) == ['$', 'E', 'l', 'e', 'm', 'e', 'n', 't', 's', '\n']:
        break

size = int(file.readline())
A = numpy.zeros(shape=(size,1))

for i in range(size):
    f = file.readline()
    f = numpy.fromstring(f, dtype=int, sep=' ')
    A[i] = f[6]

nbSsDomains = int(numpy.amax(A))
print "*************************************************************"
print "**                 INFOS PYTHON                              "
print "**  Le fichier lu est : ", filename
print "**  Le nombre de sous-domaines est :", nbSsDomains
print "**  Ce nombre est stocké dans le fichier : " "python_res.txt"
print "*************************************************************"
print
print
print

resultats = open("python_res.txt", "w")
resultats.write(filename)
resultats.write('\n')
resultats.write(str(nbSsDomains))


resultats.close()

file.close()

# Lancer le programme : mpirun -np nbSsDomains+1 exe
os.system("mpirun -np {0} exe".format(nbSsDomains+1)) 
