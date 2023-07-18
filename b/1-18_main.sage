RIFpi = RIF.pi() #4*arctan(RIF(1))
I = RIF(1)

from sage.symbolic.expression_conversions import RingConverter
SR2RIF = RingConverter(RIF)

the_path=""

# Input from the table 
load(the_path+'library/input_ternary.sage')

# All calculatory functions are here
load(the_path+'library/small_functions_ternary.sage')

# Basic potentials functions
load(the_path+'library/potentials_basic_functions.sage')

# All drawing functions are here
load(the_path+'library/drawing.sage')

import time

#==============================  Global Variables ==============================
# Computed in main()

CASE = 0
r,s = (1,1)
ttriangles = []
ttriangles_symb = ['111', '11r', '11s', '1rr', '1rs', '1ss', 'rrr', 'rrs', 'rss', 'sss']
coronas = []
d_opt= 0
pos={}
m1,mr,ms = 0,0,0
LQ = {}
z1,zr,zs = 100,100,100
V2 = [0]*18
V=lambda i,j,k: V2[pos.get((i,j,k))]

triangles_nums = {x:0 for x in ttriangles_symb}
vp_symb= ['111', '11r', '11s', '1r1', '1rr', '1rs', '1s1', '1sr', '1ss', 'r1r', 'r1s', 'rrr', 'rrs', 'rsr', 'rss', 's1s', 'srs', 'sss']
var(*['v'+x for x in vp_symb])
vp_vars = [] # used variables for vertex potentials0
bin_c = 0

#============================= 1-18 specials =======================================
var('R')
bin_d=[
ZZ[R](27*R^4 + 112*R^3 + 62*R^2 + 72*R - 29).roots(AA)[1][0],
ZZ[R](R^8 - 4590*R^6 - 82440*R^5 + 486999*R^4 - 1938708*R^3 + 2158839*R^2 - 1312200*R + 243081).roots(AA)[2][0],
ZZ[R](1024*R^3 - 692*R^2 + 448*R - 97).roots(AA)[0][0],
ZZ[R](2*R^2 - 4*R + 1).roots(AA)[0][0],
ZZ[R](944784*R^4 - 3919104*R^3 - 2191320*R^2 - 1632960*R + 757681).roots(AA)[0][0],
ZZ[R](144*R^4 + 9216*R^3 + 133224*R^2 - 127104*R + 25633).roots(AA)[2][0],
ZZ[R](4096*R^4 + 2924*R^2 - 289).roots(AA)[1][0],
ZZ[R](108*R^2 + 288*R - 97).roots(AA)[1][0],
ZZ[R](144*R^4 - 4162200*R^2 + 390625).roots(AA)[2][0]]

bin_case = {
    1:(8, 's1'),
    2:(8, 's1'),
    3:(8, 's1'),
    4:(8, 's1'),
    5:(8, 's1'),
    6:(4, 's1'),
    7:(7, 's1'),
    8:(7, 's1'),
    9:(7, 's1'),
    10:(9, 's1'),
    11:(9, 's1'),
    12:(9, 's1'),
    13:(9, 's1'),
    14:(9, 's1'),
    15:(9, 's1'),
    16:(4, 'r1'),
    17:(5, 's1'),
    18:(5, 's1'),
    19:(6, 's1')
}

def bin_corona(c):
    bin_num, rtype = bin_case[c]
    if bin_num == 4:
        return(('s', "1111"))
    elif bin_num == 5:
        return(('s', "11sss"))
    elif bin_num == 6:
        return(('s', "1ss1s"))
    elif bin_num == 7:
        return(('s', "111s"))
    elif bin_num == 8:
        return(('s', "111"))
    elif bin_num == 9:
        return(('s', "11ss"))


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

#============================== (1-18) Potentials of tight triangles: V_xyz ========================

def matrix_vector_1_18():
    '''Given 2 coronas, returns the equation coefficients matrix and the
    right-part vector for vertex potentials '''
    global vp_vars
    m,vect= [],[]
    # sum of the vertex potentials = excess in each of 10 tight triangles
    for (x,y,z) in ttriangles:
        m.append(ttriangle_to_row(x,y,z))
        vect.append(E(z+y,x+z,x+y,x,y,z))
    # sum of the vertex potentials in each corona = 0
    # 1 xorona
    for (x,corona) in coronas:
        m.append(corona_to_row(x,corona))
        vect.append(0)

        
    # 7 lines : to take the vertex potential of the 6
    # isoceles triangles as parameters: list vpparams
    # !!!! free parameters !

    QVS = VectorSpace(QQ,18)
    current_subspace = QVS.subspace(m)
    u_vecs = []
    for vec in QVS.basis():
        vec_ss =  QVS.subspace([vec])
        if current_subspace.intersection(vec_ss).dimension() == 0:
            u_vecs.append(list(vec))
            m.append(list(vec))
            current_subspace = QVS.subspace(m)

    for vec in u_vecs:
        ind = vec.index(1)
        vp_name = triangle_name((list(pos.keys())[list(pos.values()).index(ind)]))
        print(vec,vp_name)
        vect.append(var('v'+vp_name))
        vp_vars.append(var('v'+vp_name))
    
    return(Matrix(m), vector(vect))

# ============================ Main =============================

def init_bin(case):
    global CASE, bin_c, r, s, d_opt, ttriangles, coronas, pos
    CASE=case
    (r,s, disc_nums_1rs, triangles_nums, precoronas, d_opt) = input_values(case) 
    print("Ternary density: " + str(n(d_opt)))
    bin_c = bin_case[case][0]
    d_opt = RIF(bin_d[bin_c-1]*RIFpi)
    r,s = RIF(r), RIF(s)
    print ("CASE: "+str(case) + ":   r=" + str(r) + ", s=" + str(s) + ", bin density=" + str(n(d_opt)))
    bcorona = bin_corona(case)
    coronas = [(r_value(bcorona[0]), [r_value(y) for y in bcorona[1]])]
    
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



def write_RIF(x):
    if parent(x) != RIF:
        return(str(x))
        #raise NameError("write_RIF expects RIF but " + str(parent(x)) + " was given")
    return("RIF"+str(x.endpoints()))

def save_proved():
    for c in proved_cases:
        init(c)
        save_ter(c)


Rprec = RealField(200)

def write_RIF(x):
    if parent(x) != RIF:
        return(str(x))
    return ("RIF("+str(Rprec(x.lower()))+","+str(Rprec(x.upper()))+")")

def qq2str(q):
    num = q.numerator()
    den = q.denominator()
    app = round(q, ndigits=5)
    return(f'$\\frac{{{num}}}{{{den}}}\\approx {app}$')

def writing():
    name = the_path+"output/"+str(CASE)+"_constants.sage"
    name_2 = the_path+"output/"+str(CASE)+"_tables.txt"
    with open(name,'w') as f:
        f.write("CASE = "+str(CASE)+"\n")
        f.write("r,s = "+write_RIF(r)+","+write_RIF(s)+"\n")
        f.write("d_opt = "+write_RIF((d_opt))+"\n")
        f.write("m1,mr,ms = "+str(m1)+","+str(mr)+","+str(ms)+"\n")
        f.write("V2 = "+str([x for x in V2])+"\n")
        f.write("eps = "+str(eps)+"\n")
        f.write("LQ = "+str(LQ)+"\n")
        f.write("z1,zr,zs = "+str(z1)+","+str(zr)+","+str(zs)+"\n")

    with open(name_2,'w') as f:
        f.write("CASE = "+str(CASE)+"\n")
        f.write("eps,m1,mr,ms : &"+str(qq2str(eps))+"&"+str(qq2str(m1))+"&"+str(qq2str(mr))+"&"+str(qq2str(ms))+"\n")
        f.write("V111,V11r,V11s,V1r1 : &"+str(qq2str(V(I,I,I)))+"&"+str(qq2str(V(I,I,r)))+"&"+str(qq2str(V(I,I,s)))+"&"+str(qq2str(V(I,r,I)))+"\n")
        f.write("V1rr,V1rs,V1s1,V1sr : &"+str(qq2str(V(I,r,r)))+"&"+str(qq2str(V(I,r,s)))+"&"+str(qq2str(V(I,s,I)))+"&"+str(qq2str(V(I,s,r)))+"\n")
        f.write("V1ss,Vr1r,Vr1s,Vrrr : &"+str(qq2str(V(I,s,s)))+"&"+str(qq2str(V(r,I,r)))+"&"+str(qq2str(V(r,I,s)))+"&"+str(qq2str(V(r,r,r)))+"\n")
        f.write("Vrrs,Vrsr,Vrss,Vs1s : &"+str(qq2str(V(r,r,s)))+"&"+str(qq2str(V(r,s,r)))+"&"+str(qq2str(V(r,s,s)))+"&"+str(qq2str(V(s,I,s)))+"\n")
        f.write("Vsrs,Vsss : &"+str(qq2str(V(s,r,s)))+"&"+str(qq2str(V(s,s,s)))+"\n")
        f.write("z1,zr,zs : "+str(qq2str(z1))+"&"+str(qq2str(zr))+"&"+str(qq2str(zs))+"\n")	
        f.write("LQ: 11,1r,1s : "+str(qq2str(LQ[(I,I)][0]))+"&"+str(qq2str(LQ[(I,I)][1]))+"&"+str(qq2str(LQ[(I,r)][0]))+"&"+str(qq2str(LQ[(I,r)][1]))+"&"+str(qq2str(LQ[(I,s)][0]))+"&"+str(qq2str(LQ[(I,s)][1]))+"\n")
        f.write("LQ: rr,rs,ss : "+str(qq2str(LQ[(r,r)][0]))+"&"+str(qq2str(LQ[(r,r)][1]))+"&"+str(qq2str(LQ[(r,s)][0]))+"&"+str(qq2str(LQ[(r,s)][1]))+"&"+str(qq2str(LQ[(s,s)][0]))+"&"+str(qq2str(LQ[(s,s)][1]))+"\n")
    print("The constants are written in "+name)
   
