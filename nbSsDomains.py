# -*-coding:utf-8 -*

import os
import sys
import numpy

filename = raw_input("Entrer le fichier mesh dont on veut extraire le nombre de sous-domaines \n")
print "__________________________________________________"
print "Le fichier lu est : ",filename

# Test de l'existence du fichier de maillage
try:
    with open(filename): pass
except IOError:
    print "_________________________________________"
    print "Erreur : le fichier spécifié n'existe pas"
    print "_________________________________________"
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
    #f = list(f)
    #print f
    #f = [g for g in f if not g.isspace()]
    #print f
    A[i] = f[6]

nbSsDomains = int(numpy.amax(A))
print "__________________________________________________"
print "Le nombre de sous-domaines est :", nbSsDomains
print "Ce nombre est stocké dans le fichier : " "python_res.txt"
print "__________________________________________________"

resultats = open("python_res.txt", "w")
resultats.write(filename)
resultats.write('\n')
resultats.write(str(nbSsDomains))


resultats.close()

file.close()

