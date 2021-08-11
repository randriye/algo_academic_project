#!/usr/bin/env python3

"""
compute sizes of all connected components.
sort and display.
"""

from timeit import timeit
from sys import argv
import time as time
from geo.point import Point
import numpy as np
from operator import itemgetter

class UnionFind:
    """
    Structure de données pour stocker et faire les unions de cases
    cases est un dictionnaire. Les clés sont les coordonnées de cases.
    Les valeurs sont des sets de points contenus dans les cases.
    """
    def __init__(self, cases):
        """ 
        ON REPERE LES CASES PAR LEURS INDICES DANS UNE LISTE L DE CLES 
        """
        L = list(cases.keys())
        N = len(L)
        self.rank = [0 for _ in range(N)]
        self.parent = [-1 for _ in range(N)]
        self.children = [[i] for i in range(N)]
        
        # indices permet d'avoir l'indice dans L correspondant à une case
        indexes = {}
        cpt = 0
        for case in L:
            indexes[case] = cpt
            cpt += 1
        self.indices = indexes
        
        # squares permet d'accéder aux coordonnées d'une case à partir de l'indice correspondant de la liste L
        coord = {}
        for i in range(N):
            coord[i] = L[i]
        self.squares = coord
    
    def find(self, x):
        """
        Permet de trouver la racine d'un élément donné.
        On travaille sur les indices des carrés dans la liste L.
        """
        parent = self.parent[x]
        while parent >= 0:
            x = parent
            parent = self.parent[x]
        return x

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)

        if root_x == root_y:
            return
        else:
            if self.rank[root_x] >= self.rank[root_y]:
                self.parent[root_x] += self.parent[root_y]
                self.parent[root_y] = root_x
                self.rank[root_x] = max(self.rank[root_y] + 1, self.rank[root_x])
                self.children[root_x] += self.children[root_y]

            else:
                self.parent[root_y] += self.parent[root_x]
                self.parent[root_x] = root_y
                self.rank[root_y] = max(self.rank[root_x] + 1, self.rank[root_y])
                self.children[root_y] += self.children[root_x]

    def same(self, x,y):
        return self.find(x) == self.find(y)
    
    def find_children(self):
        return [(idx, -val, self.children[idx]) for idx, val in enumerate(self.parent) if val < 0]

def load_instance(filename):
    """
    loads .pts file.
    returns distance limit and points.
    """
    with open(filename, "r") as instance_file:
        lines = iter(instance_file)
        distance = float(next(lines))
        points = [Point([float(f) for f in l.split(",")]) for l in lines]

    return distance, points


def dist_sq(p1, p2):
    return (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2


def trouve_case(point, cote):
    """
    Renvoie les entiers (i,j) qui repère une case (entiers du coin supérieur droit du carré)
    """
    approx_abs = point[0] / cote
    approx_ord = point[1] / cote
    i = int(approx_abs) + 1
    j = int(approx_ord) + 1
    return (i, j)


def division(points, distance):
    """
    Renvoie un dictionnaire dont les clé sont les coordonnées des cases de diagonale distance.
    Les valeurs sont des ensembles de points contenus dans chaque case.
    """
    cases = {}
    cote = distance / np.sqrt(2)
    for k in range(len(points)):
        x, y = points[k].coordinates
        pt = (x, y)
        (i, j) = trouve_case(pt, cote)

        if (i, j) in cases:
            cases[i, j].append(pt)
        else:
            cases[i, j] = [pt]
    return cases


def max_dist_bary(l, b):
    """
    Renvoie le point de la liste l le plus éloigné du point b.
    """
    L = []
    for point in l:
        L.append( [point, dist_sq(b, point)])
    L.sort(key = itemgetter(1))
    return L[-1]


def pre_filter( key1, key2, l1, l2, distance):
    """
    Pré-filtre avec le check de cases
    Regarde si les plus éloignés des barycentres de leurs cases sont suffisamment proche
    pour que les cases soient connexes (CONDITION NECESSAIRE).
    Renvoie True si la condition n'est pas respectée (True car on est bien dans le pire cas). False si elle l'est.
    """
    G1 = barycentre(key1, distance)
    G2 = barycentre(key2, distance)
    max1 = max_dist_bary(l1, G1)
    max2 = max_dist_bary(l2, G2)
    if dist_sq(max1[0], max2[0] ) > 4 * (distance ** 2):
        return True
    else:
        return False


def is_connexe(l1, l2, distance):
    """
    Renvoie True si l1 et l2 sont connexes, False Sinon.
    En général, on travail avec des listes triées par abs/ord croissantes/décroissantes.
    """
    dist2 = distance ** 2
    if l1 == [] or l2 == []:
        return False
    # G1 = barycentre(key1, distance)
    # G2 = barycentre(key2, distance)
    for p1 in l1:
        for p2 in l2:
            if dist_sq(p1, p2) <= dist2:
                return True
    return False


def barycentre(key, distance):
    """
    renvoie le barycentre des coins d'une case de diagonale distance.
    x * cote, y * cote représente le coin superieur droit de la case
    """
    x, y = key[0], key[1]
    cote = distance / np.sqrt(2)
    c1 = [x * cote, y * cote]
    c2 = [c1[0], c1[1] - cote]
    c3 = [c1[0] - cote, c1[1] - cote]
    c4 = [c1[0] - cote, c1[1]]
    b = [0, 0]
    n = 4
    for coin in [c1, c2, c3, c4]:
        b[0], b[1] = b[0] + coin[0] / n, b[1] + coin[1] / n
    return tuple(b)
    
    
def find_couples(cases, distance):
    """
    Crée une liste de tuples. Chaque tuple est un couple de cases (créées par division)
    """
    dont_check = {}
    couples = [] # Liste des couples
    for case in cases:
        i = case[0]
        j = case[1]

        cle_ns = cases[case].copy()
        cle_sorted_ord_in = cases[case].copy()
        cle_sorted_ord_in.sort(key=itemgetter(1))
        cle_sorted_ord_de = cle_sorted_ord_in.copy()
        cle_sorted_ord_de.reverse()
        
        cle_sorted_abs_in = cases[case].copy()
        cle_sorted_abs_in.sort()
        cle_sorted_abs_de = cle_sorted_abs_in.copy()
        cle_sorted_abs_de.reverse() 
        """ CAS J """
        if (i - 2, j) not in dont_check and (i - 2, j) in cases:
            to_check = cases[i - 2, j]
            worst = pre_filter((i, j), (i - 2, j), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(reverse = True)
                im2j = is_connexe(to_check, cle_sorted_abs_in, distance)

                if im2j:
                    couples.append(((i, j),(i - 2, j)))

        #couple_treated.add( ((i,j), (i - 2 *cote, j)) )
                if im2j:
                    couples.append(((i,j),(i - 2, j)))
            #couples[cle].add((i - 2 * cote, j))


        if (i - 1, j) not in dont_check and (i - 1, j) in cases:
            to_check = cases[i - 1, j]
            worst = pre_filter((i, j), (i - 1, j), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(reverse = True)
                im1j = is_connexe(to_check, cle_sorted_abs_in, distance)
                if im1j:
                    couples.append(((i, j), (i - 1, j)))

        if (i + 1, j) not in dont_check and (i + 1, j)  in cases:
            to_check = cases[i + 1, j]
            worst = pre_filter((i, j), (i + 1, j), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort()
                ip1j = is_connexe(to_check, cle_sorted_abs_de, distance)
                if ip1j:
                    couples.append( ((i, j), (i + 1, j)))

        if (i + 2 , j) not in dont_check and (i + 2 , j) in cases:
            to_check = list(cases[i + 2, j])
            worst = pre_filter( (i,j), (i + 2, j), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort()
                ip2j = is_connexe( to_check, cle_sorted_abs_de, distance)
                if ip2j:
                    couples.append( ( (i,j), (i + 2 , j)))
        
        """ CAS J - 1 """
        if (i - 2 , j - 1) not in dont_check and (i - 2 , j - 1) in cases: 
            to_check = cases[i - 2 , j - 1 ]
            worst = pre_filter( (i,j), (i - 2, j - 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(reverse=True)
                im2jm1 = is_connexe(to_check, cle_sorted_abs_in, distance)
                if im2jm1:
                    couples.append( ((i,j) ,(i - 2 , j - 1) ))
        
        if (i - 1 , j - 1 ) not in dont_check and (i - 1 , j - 1 ) in cases:
            to_check = cases[i - 1 , j - 1 ]
            worst = pre_filter((i,j), (i - 1, j - 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                im1jm1 = is_connexe(to_check, cle_sorted_abs_de, distance)
                if im1jm1:
                    couples.append( ((i,j), (i - 1 , j - 1)))
        
        if (i, j - 1 ) not in dont_check and (i, j - 1 ) in cases:
            to_check = cases[i, j - 1 ]
            
            worst = pre_filter((i,j), (i, j - 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                ijm1 = is_connexe( to_check, cle_sorted_ord_in, distance)
                if ijm1:
                    couples.append( ( (i,j), (i, j - 1 )))

        if (i + 1 , j - 1 ) not in dont_check and (i + 1 , j - 1 ) in cases:
            to_check = cases[i + 1 , j - 1 ]
            worst = pre_filter((i,j), (i + 1, j - 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                ip1jm1 = is_connexe( to_check, cle_sorted_ord_de, distance)
                if ip1jm1:
                    couples.append( ((i,j), (i + 1 , j - 1 )))

        if (i + 2, j - 1 ) not in dont_check and (i + 2, j - 1 )  in cases:
            to_check = cases[i + 2, j - 1 ]
            worst = pre_filter((i,j), (i + 2, j - 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort()
                ip2jm1 = is_connexe(to_check, cle_sorted_abs_de, distance)
                if ip2jm1:
                    couples.append( ((i,j), (i + 2, j - 1 )))

        """ CAS J - 2"""
        if (i - 1 , j - 2 ) not in dont_check  and (i - 1 , j - 2 ) in cases:
            #if ( (i - 1 * cote, j - 2 * cote), (i,j)) not in couple_treated:
            to_check = cases[i - 1 , j - 2 ]
            worst = pre_filter((i,j), (i - 1, j - 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                im1jm2 = is_connexe(to_check, cle_sorted_ord_in, distance)
            #couple_treated.add( ((i,j), (i - 1 * cote, j - 2 * cote)))
                if im1jm2:
                    couples.append( ((i,j), (i - 1 , j - 2 )))

        if (i, j - 2 ) not in dont_check and (i, j - 2 ) in cases:
            #if ( (i, j - 2 * cote), (i,j)) not in couple_treated:
            to_check = cases[i, j - 2 ]
            worst = pre_filter((i,j), (i, j- 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                ijm2 = is_connexe( to_check, cle_sorted_ord_in, distance)
            #couple_treated.add( ((i,j), ( i, j - 2 * cote)) )
                if ijm2:
                    couples.append( ( (i,j), (i, j - 2 )))

        if (i + 1 , j - 2 ) not in dont_check and (i + 1 , j - 2 ) in cases:
            #if( (i + 1 * cote, j - 2 * cote), (i,j)) not in couple_treated:
            to_check = cases[i + 1 , j - 2 ]
            worst = pre_filter((i,j), (i + 1, j - 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1), reverse=True)
                ip1jm2 = is_connexe(to_check, cle_sorted_ord_in, distance)
            #couple_treated.add( ((i,j), (i + 1 * cote, j - 2 * cote)))
                if ip1jm2:
                    couples.append( ((i,j), (i + 1 , j - 2 )))

        """CAS J + 1"""
        if (i - 2 , j + 1 ) not in dont_check and (i - 2 , j + 1 ) in cases:
            #if ( (i - 2 * cote, j + 1 * cote), (i,j)) not in couple_treated:
            to_check = cases[i - 2 , j + 1 ]
            worst = pre_filter((i,j), (i - 2, j + 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(reverse=True)
                im2jp1 = is_connexe( to_check, cle_sorted_abs_in, distance)
            #couple_treated.add( ((i,j), (i - 2 * cote, j + 1 * cote)))
                if im2jp1:
                    couples.append( ((i,j), (i - 2 , j + 1 )))

        if (i - 1 , j + 1) not in dont_check and (i - 1 , j + 1) in cases:
            #if ( (i - 1 * cote, j + 1 * cote), (i,j)) not in couple_treated:
            #print("on est dedans", composantes[i,j])
            to_check = cases[i - 1 , j + 1 ]
            worst = pre_filter( (i,j), (i - 2, j + 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                im1jp1 = is_connexe( to_check, cle_sorted_ord_de, distance)
            #couple_treated.add( ((i,j), (i - 1 * cote, j + 1 * cote)))
                if im1jp1:
                    couples.append( ((i,j), (i - 1 , j + 1)))

        if (i, j + 1) not in dont_check and (i, j + 1) in cases:
            to_check = cases[i, j + 1]
            worst = pre_filter( (i,j), (i - 2, j + 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort( key = itemgetter(1))
                ijp1 = is_connexe ( to_check, cle_sorted_ord_de, distance)
                if ijp1:
                    couples.append( ((i,j), (i, j + 1)))



        if (i + 1 , j + 1 ) not in dont_check and (i + 1 , j + 1 ) in cases:
            #if ( (i + 1 * cote, j + 1 * cote), (i,j)) not in couple_treated:
            to_check = cases[i + 1 , j + 1 ]
            worst = pre_filter((i,j), (i + 1, j + 1), cle_ns,to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                ip1jp1 = is_connexe( to_check, cle_sorted_ord_de, distance)
            #couple_treated.add( ((i,j), (i + 1 * cote, j + 1 * cote)))
                if ip1jp1:
                    couples.append( ((i,j), (i + 1 , j + 1 )))

        if (i + 2 , j + 1 ) not in dont_check and (i + 2 , j + 1 ) in cases:
            #if ( (i + 2 * cote, j + 1 * cote), (i,j)) not in couple_treated:
            to_check = cases[i + 2 , j + 1 ]
            worst = pre_filter((i,j), (i + 2, j + 1), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                ip2jp1 = is_connexe( to_check, cle_sorted_ord_de, distance)
            #couple_treated.add ( ((i,j), (i + 2 * cote, j + 1 * cote)))
                if ip2jp1:
                    couples.append( ((i,j), (i + 2 , j + 1 )))
        
        """ CAS J + 2 """

        if (i - 1 , j + 2 ) not in dont_check and (i - 1 , j + 2 ) in cases:
            #if ( (i - 1 * cote, j + 2 * cote), (i,j)) not in couple_treated:
            to_check = cases[i - 1 , j + 2 ]
            worst = pre_filter( (i,j), (i - 1, j + 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                im1jp2 = is_connexe( to_check, cle_sorted_ord_de, distance)
            #couple_treated.add ( ((i,j), (i - 1 * cote, j + 2 * cote)))
                if im1jp2:
                    couples.append( ((i,j), (i - 1 , j + 2 )))

        if (i, j + 2 ) not in dont_check and (i, j + 2 ) in cases:
            #if ( (i, j - 2 * cote), (i,j)) not in couple_treated:
            to_check = cases[i, j + 2 ]
            worst = pre_filter((i,j), (i,j + 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                ijm2 = is_connexe( to_check, cle_sorted_ord_de, distance)
                if ijm2:
                    couples.append( ((i,j), (i, j + 2)))

        if (i + 1 , j + 2 ) not in dont_check and (i + 1 , j + 2 ) in cases:
            to_check = cases[i + 1 , j + 2 ]
            worst = pre_filter( (i,j), (i + 1, j + 2), cle_ns, to_check, distance)
            if worst == False:
                to_check.sort(key=itemgetter(1))
                ip1jp2 = is_connexe( to_check, cle_sorted_ord_de, distance)
                if ip1jp2:
                    couples.append( ((i,j), (i + 1 , j + 2 )))
        
        dont_check[i,j] = True
    return couples

def fusion(couples, cases):
    """
    Prend une pile/liste de couples de cases connexes 2 à 2, et les unifie dans un UnionFind.
    S'arrếte quand la pile est vide.
    """
    composantes = UnionFind(cases)
    while couples != [] :
        couple = couples.pop()
        composantes.union(composantes.indices[couple[0]], composantes.indices[couple[1]])
    return composantes

def compute_length(composantes, cases):
    """
    Récupère un UnionFind et calcule les tailles des composantes.
    Utilise la méthode find_children de la classe UnionFind.
    On prend tous les children d'une composante et on prend leur taille pour les sommer.
    """
    L = composantes.find_children()
    lengths = [] 
    for elem in L:
        children = elem[2]
        length = 0
        for child in children:
            length += len( cases[composantes.squares[child]])
        lengths.append(length)
    return lengths


def print_components_sizes(distance, points):
    """
    Affiche la taille des composantes de la liste points avec la contrainte distance, par ordre décroissant.
    """
    if distance <= 0:
        # CAS TOUS LES POINTS SONT SEULS
        print([1 for _ in range(len(points))])
    elif distance >= np.sqrt(2):
        # CAS TOUS LES POINTS SONT ENSEMBLES
        print([len(points)])
    else:
        # CAS GENERAL

        #dd = time.time()
        cases = division(points, distance)
        #fd = time.time()
        #print("division =", fd - dd, "s")
        #dfc = time.time()
        couples = find_couples(cases, distance)
        #ffc = time.time()
        #print( "find couples", ffc - dfc, "s")
        #df = time.time()
        composantes = fusion(couples, cases)
        #ff = time.time()
        #print("fusion", ff - df, "s")
        #dl = time.time()
        lengths = compute_length(composantes, cases)
        lengths.sort(reverse = True)
        #fl = time.time()
        #print("length = ", fl - dl, "s")
        print(lengths)


def main():
    """
    ne pas modifier: on charge une instance et on affiche les tailles
    """
    for instance in argv[1:]:
        distance, points = load_instance(instance)
        print_components_sizes(distance, points)


main()
