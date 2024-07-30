
def MakeChampABSC_CURV(mesh, grma):
    '''
    Return a field of curvilinear coordinates for a specific group of elements.
    Input:
        mesh: mesh to extract field from.
        grma: name of the elements' group.
    Output:
        chacno: nodal field of curvilinear coordinates.
    '''
    chabsc = CREA_CHAMP(identifier='9:1',
                        AFFE_SP=_F(CARA_ELEM=elprop),
                        INFO=1,
                        MAILLAGE=mesh,
                        NOM_CHAM='ABSC_CURV',
                        OPERATION='EXTR',
                        TYPE_CHAM='CART_ABSC_R')

    tabsc = CREA_TABLE(identifier='10:1',
                      RESU=_F(CHAM_GD=chabsc,
                              GROUP_MA=grma,
                              TOUT_CMP='OUI')).EXTR_TABLE()

    chgnoeu = CREA_CHAMP(identifier='11:1',
                         MAILLAGE=mesh,
                         NOM_CHAM='GEOMETRIE',
                         OPERATION='EXTR',
                         TYPE_CHAM='NOEU_GEOM_R')

    chgelno = CREA_CHAMP(identifier='12:1',
                         CHAM_GD=chgnoeu,
                         INFO=1,
                         MODELE=model,
                         OPERATION='DISC',
                         TYPE_CHAM='ELNO_GEOM_R')

    telno = CREA_TABLE(identifier='13:1',
                      RESU=_F(CHAM_GD=chgelno,
                              GROUP_MA=grma,
                              TOUT_CMP='OUI')).EXTR_TABLE()

    valac = num.vstack([list(tabsc.ABSC1),
                        list(tabsc.ABSC2),
                        list(tabsc.ABSC3)]).flatten('F')
    
    noeac = list(telno.NOEUD)

    vord, pord = num.unique(valac, return_index=True)
    vord = vord.tolist()
    nord = num.array(noeac)[pord].tolist()

    tacno = CREA_TABLE(
        LISTE=(
            _F(LISTE_K=nord, TYPE_K='K8', PARA='NOEUD'),
            _F(LISTE_R=vord, PARA='X1'),
        ),
    )

    chacno = CREA_CHAMP(
        TYPE_CHAM='NOEU_NEUT_R', OPERATION='EXTR', TABLE=tacno, MAILLAGE=mesh
    )

    return chacno

def dfs_elem(C, P, N, eo, n, em):
    '''
    Separate linear elements in members.
    This is a recursive function.
    Input:
        C: connectivity matrix (elem -> [nodes]).
        P: pertinence matrix (node -> [elem]).
        N: nodes' coordinates.
        eo: starting element.
        n: member family at present level.
        em: element-member matrix.
    '''
    # Percorrer os nós do elemento.
    for i in C[eo]:
        # Se for uma conexão simples,
        if len(P[i]) <= 2:
              # buscar nos elementos que contém o nó.
              for ex in P[i]:
                  # Se não tiver sido marcado, for diferente do original e for 
                  # colinear.
                  vo = N[C[eo][0]] - N[C[eo][1]]
                  vx = N[C[ex][0]] - N[C[ex][1]]
                  if ex != eo and em[ex] == -1 and \
                                   num.linalg.norm(num.cross(vo, vx)) < 1e-10:
                      # Marcar e passar a nova busca com base neste.
                      em[ex] = n
                      em = dfs_elem(C, P, N, ex, n, em)
    return (em)
    
def node_elem(C, N):
    '''
    Make a pertinence matrix, relating a node to the elements if pertains to.
    Input:
        C: connectivity matrix.
        N: nodes.
    Output:
        P: pertinence
    '''
    P = dict(zip( N.keys(), [[]]*len(N)))
    for e,n in C.items():
        for i in C[e]:
            if P[i] == []:
                P[i] = [e]
            else:
                P[i].append(e)
    return(P)

def absc_curv(C, P, V, norig, grma):
    '''
    Calculates the curvilinear coordinates of a specific element family.
    Input:
        C: connectivity matrix.
        P: pertinence matrix.
        V: vertexes (nodes) coordinates.
        norig: first node.
        grma: element group to sweep.
    Output:
        AC: curvilinear coordinates.
        ndone: list of nodes.
    '''
    etodo = list(grma)
    edone = []
    ntodo = list(P.keys())
    ndone = []
    
    # Initial node for summing the curvilinear coordinate.
    n = norig
    
    # Sweeping through the connections.
    while len(etodo) > 0:
        # Finding the element to sweep.
        for e in P[n]:
            if e in etodo:
                break
        
        blk = []
        if n == C[e][0]: # The node is the first one in the element.
            blk = [C[e][0], *C[e][2:], C[e][1]]
        elif n == C[e][1]: # The node is the last one in the element.
            blk = [C[e][0], *C[e][2:], C[e][1]]
            blk.reverse()
        else:
            print('Failed:{0}',e)
    
        # Next point.
        n = blk[-1]

        # Erasing duplicates.        
        for nx in blk:
            if ndone.count(nx) != 0:
                blk.pop(blk.index(nx))
        
        # Appending the block of ordered nodes.
        ndone = ndone + blk
    
        # Marking the nodes as 'done'.
        for nx in blk:
            if nx in ntodo:
                ntodo.pop(ntodo.index(nx))

        # Marking the element as 'done'.
        edone.append(etodo.pop(etodo.index(e)))
        
    # Accumulating the curvilinear coordinates in the order.
    AC = {}
    AC[norig] = 0
    blk = len(ndone)
    for i in range(1,len(ndone)):
        AC[ndone[i]] = AC[ndone[i - 1]] \
            + num.linalg.norm(V[ndone[i]] - V[ndone[i - 1]])
        
    return(AC, ndone)
