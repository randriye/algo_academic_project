Academic Project 

Objective : Find the sizes of connected components in a scatter graph

We are given a list of points in [0,1]x[0,1]. The objective is to return a list of the sizes of the connected components of those points considering a certain threshold distance we are going to call "s".

Team : Group of 2
Time Spent : ~40h
Technologies used:

- Python
- Libraries : geo, numpy, matplotlib
- Structures : python-list, dictionnaries, union-find
- Mathematics : Euclidian Distance, Triangular inequality, Barycenters

Details of the Method:
- Divide the square [0,1]x[0,1] in smaller squares of side s / sqrt(2) => every non-empty square represents now a connected components. The components are stocked in a dictionnary : keys are coordinates of the squares, values are lists of points of each component.
- For each component, we check the 20 cases around it (all cases except the corners). To accelerate, we sort by abscissas or ordinates, and by ascending or descending order. The triangular inequality  and the barycenter are used to check the case if the components are connected or not. When a case is treated, it is not checked by the other cases.
- To stock the connection between cases, we use union-find structure, which has an average cost of O(1) in Research and in Union. The worst case is in O(log n). Each case is indexed by a unique number.
- Browse all big connected components (after the union) to check its size (the number of points)
- Print the the sorted list
