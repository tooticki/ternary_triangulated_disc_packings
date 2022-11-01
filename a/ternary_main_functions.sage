import time

the_path=""

# Input from the table 
load(the_path+'library/input_ternary.sage')

# All calculatory functions are here
load(the_path+'library/small_functions_ternary.sage')

# Basic potentials functions
load(the_path+'library/potentials_basic_functions.sage')

# All drawing functions are here
load(the_path+'library/drawing.sage')

#==============================  Global Variables ==============================

RIFpi = 4*arctan(RIF(1)) # Correct value of RIF(pi)
I = RIF(1)
ttriangles_symb = ['111', '11r', '11s', '1rr', '1rs', '1ss', 'rrr', 'rrs', 'rss', 'sss']

# Computed in main
CASE = 0
r,s = (1,1)
ttriangles = []
coronas = []
d_opt= 0
pos={}
V2 = [0]*18
m1,mr,ms = 0,0,0
LQ = {}
z1,zr,zs = 100,100,100
#V(i,j,k) = vertex potential V_{ijk}
V=lambda i,j,k: V2[pos.get((i,j,k))]



#============================== Potentials of tight triangles: V_xyz ==============================

def ttriangle_to_row(x,y,z):
    '''Given the radii of circles of a tight triangle, returns the raw in
    the vertex-potential matrix corresponding to this triangle'''
    row = [0]*18
    row[pos[(x,y,z)]]+=1
    row[pos[(y,z,x)]]+=1
    row[pos[(z,x,y)]]+=1
    return row

def corona_to_row(x, corona):
    '''Given the radius of the circle x and its corona, returns the raw in the
    vertex-potential matrix corresponding to this corona  '''
    row = [0]*18
    n = len(corona)
    for i in range(n):
        y = corona[i-1]
        z = corona[i]
        row[pos[(y,x,z)]]+=1
    return row

def isvertice_to_row(v):
    '''Returns a row of the matrix corresponding to the isosceles triangle
    with radii y,x,y where (x,y)=v'''
    (x,y) = v
    row = [0]*18
    row[pos[(y,x,y)]]+=1
    return row

var('vr1r','vs1s','v1r1','vsrs','v1s1','vrsr') # isvertices trianlges' potentials

# ============================== Main =============================

def matrix_vector(coronas):
    '''Given 2 coronas, returns the equation coefficients matrix and the
    right-part vector for vertex potentials '''
    matrix,vect= [],[]
    # sum of the vertex potentials = excess in each of 10 tight triangles
    for (x,y,z) in ttriangles:
        matrix.append(ttriangle_to_row(x,y,z))
        vect.append(E(z+y,x+z,x+y,x,y,z))
    # sum of the vertex potentials in each corona = 0
    # (there are always <= 4 different coronas) normally, 2
    for (x,corona) in coronas:
        matrix.append(corona_to_row(x,corona))
        vect.append(0)
    # 6 lines (or <= 6) : to take the vertex potential of the 6
    # isoceles triangles as parameters: list vpparams
    # !!!! free parameters !
    isvertices = [(x,y) for x in [I,r,s] for y in [I,r,s] if x!=y]
    vpparams = {(I,r):vr1r,
                (I,s):vs1s,
                (r,I):v1r1,
                (r,s):vsrs,
                (s,I):v1s1,
                (s,r):vrsr}
    for v in isvertices:
        matrix.append(isvertice_to_row(v))
        vect.append(vpparams[v])
    return(Matrix(matrix), vector(vect))



def init(case):
    global CASE, r, s, d_opt, ttriangles, coronas, pos
    CASE=case    
    (r,s, disc_nums_1rs, triangles_nums, precoronas, d_opt) = input_values(case)

    r,s = RIF(r), RIF(s)
    print ("CASE "+str(case) + ":   r=" + str(r) + ", s=" + str(s) + ", density=" + str(n(d_opt)))

    coronas = [(s, [r_value(y) for y in precoronas[0]]), (r, [r_value(y) for y in precoronas[1]])]
    ttriangles = [(I,I,I),(I,I,r),(I,I,s),(I,r,r),(I,r,s),(I,s,s),(r,r,r),(r,r,s),(r,s,s),(s,s,s)]

    pos = {
        (I,I,I):0,
        (r,r,r):1,
        (s,s,s):2,
        (s,I,s):3,
        (I,s,s):4,
        (s,s,I):4,
        (r,I,r):5,
        (I,r,r):6,
        (r,r,I):6,
        (s,r,s):7,
        (r,s,s):8,
        (s,s,r):8,
        (I,s,I):9,
        (I,I,s):10,
        (s,I,I):10,
        (I,r,I):11,
        (I,I,r):12,
        (r,I,I):12,
        (r,s,r):13,
        (r,r,s):14,
        (s,r,r):14,
        (I,s,r):15,
        (r,s,I):15,
        (r,I,s):16,
        (s,I,r):16,
        (I,r,s):17,
        (s,r,I):17}

# ============================== Resluts =============================
proved_cases_a = [53,54,55,56,66,76,77,79,93,108,115,116,118,129,131, 146]
proved_cases_b = list(range(1,16))
not_saturated = [24,28,29,30,31,32,33,37,38,39,40,41,42,43,44]
counter_ex =  [19,20,25,47,51,60,63,64,70,73,80,92,95,97,98,99,100,104,110,111,117,119,126,132,133,135,136,137,138,139,141,142,151,152,154,159,161,162,163,164]
two_coronas = [16,17,18,25,31,36,39,43,49,52,57,58,65,73,78,84,90,99,106,110,111,114,117,120,137,142,148]+list(range(153,160,1))+[162,164]
empty_polyhedra =[21,22,23,26,27,33,34,35,46,48,50,59,61,67,68,69,71,72,74,81,82,83,85,86,87,88,89,91,94,96,101,102,103,105,107,109,112,113,121,122,123,124,125,127,128,130,134,140,143,145,147,149,150,160]
mysterious = [45,62,75,144]

