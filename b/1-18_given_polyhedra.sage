the_path = ""

load(the_path+"1-18_main.sage")

import time

#============================== Find epsilon > 0 ==============================

def find_epsilon_rec(a, b, lb1, lbr, lbs):
    ''' Recurrence : true on a and false on b'''
    if b-a < 0.0000000001:
        if a>0:
            return a
        else:
            NameError("No positive epsilon works with these Vs and ms")
    c = (a+b)/2
    if test_eps(lb1,lbr,lbs,c):
        return find_epsilon_rec(c,b,lb1,lbr,lbs)
    else:
        return find_epsilon_rec(a,c,lb1,lbr,lbs)
    

def find_epsilon(lb1, lbr, lbs, epsmax):
    ''' Looks for eps using dichotomic search...'''
    if test_eps(lb1, lbr, lbs, epsmax):  
        return epsmax 
    if not test_eps(lb1, lbr, lbs, QQ(0)):
        raise NameError("No positive epsilon works with these Vs and ms")
    return find_epsilon_rec(QQ(0), epsmax, lb1, lbr, lbs)


def test_eps(mm1,mmr,mms,eps):
    ''' Tests if for the given eps, M1,Mr,Ms satisfy min de >= max du'''
    var('a','b','c')
    for (Ra,Rb,Rc) in ttriangles:
        for x in [a,b,c]:
            du=RIF(dU(Ra,Rb,Rc,eps,x).subs({M1:mm1,Mr:mmr,Ms:mms}))
            de=RIF(dE(Ra,Rb,Rc,eps,x).subs({M1:mm1,Mr:mmr,Ms:mms}))
            if du.upper() > de.lower():
                return False
    return True

#============================== m-polyhedron ==============================

def choose_point(P):
    ''' A point between the vertex minimizing m1+mr+ms and the center'''
    weight = QQ(.999) #  the weight of the found point
    v_minm1 = min(P.vertices(), key = lambda v: v[-3])
    v_minmr = min(P.vertices(), key = lambda v: v[-2])
    v_minms = min(P.vertices(), key = lambda v: v[-1])
    v_minzs =  max(P.vertices(), key = lambda v: v[-5])
    print([n(x) for x in v_minzs]) 
    v_minm = [(x1+xr+4*xs)/6 for (x1,xr,xs) in zip(v_minm1, v_minmr, v_minms)] #1 1 4 normally
    center = P.center()
    the_point = [weight*v + (1-weight)*c for (v,c) in zip(v_minm, center)]
    return the_point

def verify_vm(q):
    '''since we reduced the number of inequalities to build Vmp1,Vmpr,Vmps, we
    still need to make sure that the choose m's and vijk satisfy ALL the
    m-inequalities. That is why we use this function'''
    print("    Checking m"+r_name(q)+"...")
    mq = mq_value(q)
    timer = time.time()
    for triangles in corona_generator(q): #only for 10
        sum_v = sum(V(*t) for t in triangles)
        
        sum_angles = sum(tight_angle(t) for t in triangles)
        
        diff = 2*RIFpi - sum_angles
        
        if not diff.contains_zero() and not mq >= - sum_v/abs(diff):
            print(" !!! M"+str(r_name(q))+" = " + str(n(mq)) +  " < -sumv/diff = "+ str(n(- RIF(sum_v/abs(diff))))+" for "+str([triangle_name(t) for t in triangles]))
            return False
    print("     Checking m"+r_name(q)+ " took " + str(time.time()-timer))
    return True

#============================== Edge potentials ==============================

def find_max_l(R,x,y):
    distxy = lambda z,u: dist(z,x+u,y+u,u,y,x,radius(z,x+u,y+u,u,y,x))
    EminusU = lambda t: E(t,x+R,y+R,R,y,x) - U(t,x+R,y+R,R,y,x)
    a,b = x+y+QQ(.00000001), sqrt((R+x)**2-R**2) + sqrt((R+y)**2-R**2)
    if EminusU(a) <= 0:
        NameError (CASE + ":  Epsilon seems to be < 1e-7")
    if EminusU(b) >= 0:
        return b
    lxy = sign_change_point(EminusU, a, b, QQ(1e-4)) - QQ(.00001)
   
    return lxy

def find_q(l,R,x,y):
    a,b = max(l, abs(x-y)),  min(x+y+2*R-QQ(0.000001), sqrt((R+x)**2-R**2) + sqrt((R+y)**2-R**2)+s/QQ(2))
   
    ERxy = lambda t: E(t,x+R,y+R,R,y,x)
    E0 = sign_change_point(ERxy, a, b, QQ(1e-4)) + QQ(.00001)
    Ue0 =U(E0,x+R,y+R,R,y,x)
    re0 = radius(E0,x+R,y+R,R,y,x)
    de0 = dist(E0,x+R,y+R,R,y,x,re0)
    qxy = (-Ue0-QQ(.001))/de0

    Eb = ERxy(b)
    Ub = U(b,x+R,y+R,R,y,x)
    reb = radius(b,x+R,y+R,R,y,x)
    deb = dist(b,x+R,y+R,R,y,x,reb)
    qb =  (Eb-Ub-QQ(.001))/deb
    return max(qxy, qb)

def lq(x,y,EPS=QQ(1e-1)):
    ''' Ue(T) = qxy*de  if  xy > lxy else 0 '''
    l = QQ(min([find_max_l(R,x,y) for R in [I,r,s]]).lower()) 
    q = QQ(max([find_q(l,R,x,y) for R in [I,r,s]]).lower())    
    return (l,q)


#============================== Z1,Zr,Zs ==============================

#------------------------------ Checking that mq >= U(T)/v(T) ------------------------------

def mq_ineq(q):
    '''Return True iff for all tight T, mq >= U(t)/v(T) where v(T) is the
    angle in q-circle center : this inequality is essential for capRIFping'''
    print(q)
    mq = mq_value(q)
    for t in tripl(q):
        if mq < -V(*t)/tight_angle(t):
            print("mq <  U(T)/v(T) for T =  " + str(triangle_name(t)))
            return False
    return(True)

zq = lambda q : QQ( max([  (-2 * RIFpi * V(*t)/tight_angle(t)).upper() for t in tripl(q) ]))

#============================== Testing on all triangles ==============================

def all_triangles():
    '''Returns a list of alll triangles we need to consider'''
    # an FM-triangle should at least contain half an s-disc -> lower bound on 16 times its squared area
    Smin=(16*((RIFpi/2)*s^2)^2).lower()

    T = []
    for (x,y,z) in ttriangles:
        if Q_inequality(x,y,z): # With it, we are allowed to consider only those with 1 contact
            T += [(y+z,RIF(x+z,x+z+2*s),RIF(x+y,x+y+2*s),x,y,z)]
            T += [(RIF(y+z,y+z+2*s),x+z,RIF(x+y,x+y+2*s),x,y,z)]*(r_name(x)!=r_name(y))
            T += [(RIF(y+z,y+z+2*s),RIF(x+z,x+z+2*s),x+y,x,y,z)]*(r_name(x)!=r_name(z) and r_name(y)!=r_name(z))
        else:
            print("    Q-inequality doesn't hold on "+str((x,y,z)))
            T += [(RIF(y+z,y+z+2*s),RIF(x+z,x+z+2*s),RIF(x+y,x+y+2*s),x,y,z)]
    return T

nfeasable_num = {t:0 for t in ttriangles_symb}
epstight_num = {t:0 for t in ttriangles_symb}
good_num = {t:0 for t in ttriangles_symb}

epstight = lambda a,b,c,ra,rb,rc : a<=rb+rc+eps and b<=ra+rc+eps and c<=ra+rb+eps

# not feasible because too flat: triangle should at least contain half an s-disc
nfeasable = lambda a,b,c : Ss(a,b,c) < QQ((QQ(16)*((RIFpi/QQ(2))*s^2)^2).lower())

def test_triangle(a,b,c,ra,rb,rc):
    global nfeasable_num, epstight_num, good_num
    t_name = triangle_name((ra,rb,rc))
    if epstight(a,b,c,ra,rb,rc): 
        epstight_num[t_name]+=1
        return 0
    if nfeasable(a,b,c):  
        nfeasable_num[t_name]+=1
        return 0
    S=Ss(a,b,c)
    A=Aa(a,b,c,ra,rb,rc)
    B=Bb(a,b,c,ra,rb,rc)
    C=Cc(a,b,c,ra,rb,rc)
    D=Dd(a,b,c,ra,rb,rc)
    if S >= 0 and (not A.contains_zero() or (not B.contains_zero() and not D.contains_zero())):
        R=radiusABCD(a,b,c,ra,rb,rc,A,B,C,D)
        if R >= s:    # not feasible because not saturated
            nfeasable_num[t_name]+=1
            return 0
        ee=E(a,b,c,ra,rb,rc)
        uu=U2withR(a,b,c,ra,rb,rc,R)
        if ee > uu:   # good triangle 
            good_num[t_name]+=1
            return 0
        if ee < uu and R < s:   # bad triangle (shall not happen)
            save_triangle(a,b,c,ra,rb,rc,"__"+str(CASE)+"_BAD_")
            infos = str([x.endpoints() for x in [a,b,c]])+ " R = " +str(R.endpoints())
            raise NameError("(case "+str(CASE)+") "+t_name+" E<U: "+str(n(ee))+" < "+str(uu)+" with "+infos)

    # The actual precision does not allow to conclude -> refine
    # taking into account that some values out of a,b,c are not really
    # intervals: we should not devide them 
    EPS = 1e-9 #heuristics
    
    glue_ends = lambda point: [point] if point.diameter()<=EPS else [RIF(z) for z in list(RIF(point).bisection())]
    new_triples = cartesian_product([glue_ends(x) for x in [a,b,c]])
    
    if len(new_triples) == 1: #all of them are points: we are blocked here
        save_triangle(a,b,c,ra,rb,rc,"_blocked")
        raise NameError("(case "+str(CASE)+") "+"We're blocked on: " +t_name+" "+str((a,b,c))) 

    for (a2,b2,c2) in new_triples:
        test_triangle(a2,b2,c2,ra,rb,rc)
    return 0
    
    
def test_all():
    T = all_triangles()
    for (a,b,c,ra,rb,rc) in T:
        print(triangle_name((ra,rb,rc))+str(": testing..."))
        test_triangle(a,b,c,ra,rb,rc)
    print("---- All triangles test OK ----")
    print("Non Feasable: "+str(nfeasable_num))
    print("Eps-tight: "+str(epstight_num))
    print("Good: "+str(good_num))

def test_one(r1,r2,r3):
    tname =  triangle_name((r1,r2,r3))
    T = all_triangles()
    for (a,b,c,ra,rb,rc) in T:
        if triangle_name((ra,rb,rc)) == tname:
            test_triangle(a,b,c,ra,rb,rc)
    print("---- " + tname + " test OK ----")

#============================== Find epsilon > 0 ==============================

def find_epsilon_rec(a, b, lb1, lbr, lbs):
    ''' Recurrence : true on a and false on b'''
    if b-a < 0.0000000001:
        if a>0:
            return a
        else:
            NameError("No positive epsilon works with these Vs and ms")
    c = (a+b)/2
    if test_eps(lb1,lbr,lbs,c):
        return find_epsilon_rec(c,b,lb1,lbr,lbs)
    else:
        return find_epsilon_rec(a,c,lb1,lbr,lbs)
    

def find_epsilon(lb1, lbr, lbs, epsmax):
    ''' Looks for eps using dichotomic search...'''
    if test_eps(lb1, lbr, lbs, epsmax):  
        return epsmax # too strange
    if not test_eps(lb1, lbr, lbs, QQ(0)):
        raise NameError("No positive epsilon works with these Vs and ms")
    return find_epsilon_rec(QQ(0), epsmax, lb1, lbr, lbs)


def test_eps(mm1,mmr,mms,eps):
    ''' Tests if for the given eps, M1,Mr,Ms satisfy min de >= max du'''
    var('a','b','c')
    for (Ra,Rb,Rc) in ttriangles:
        for x in [a,b,c]:
            du=RIF(dU(Ra,Rb,Rc,eps,x).subs({M1:mm1,Mr:mmr,Ms:mms}))
            de=RIF(dE(Ra,Rb,Rc,eps,x).subs({M1:mm1,Mr:mmr,Ms:mms}))
            if du.upper() > de.lower():
                return False
    return True

#============================== Main =============================

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



        

def get_given_polyhedron(q, case):
    
    Pq = load(the_path+"output/"+str(case)+"_P_"+q+".sobj") #TODO upper
    print("  Found P"+q+": "+str(Pq))
    return Pq



def main(case):
    global V2, m1, mr, ms, eps, z1, zr, zs, LQ
    init_bin(case)
    
    # ------------------------- Computing Vxyz  -------------------------
    M,v=matrix_vector_1_18()
    V2=M.inverse()*v

    #  ------------------------- Polyhedra  -------------------------
    Pvmeps = get_given_polyhedron("e1rs",case)

    #  ------------------------- Choose a point Vm -------------------------
    vm_point = choose_point(Pvmeps) # 7 free variables
    m1,mr,ms = vm_point[-3:]
    
    v_dict = {v:val for v,val in zip(vp_vars, vm_point[:-3])}
    V2 = V2.subs(v_dict) # get rid of variables
    V2 = [QQ(RIF(x).center()) for x in V2]
    print(list(map((lambda x: (x,n3(v_dict[x]))),v_dict)))
     
    print(" m1="+str(n3(m1))+", mr="+str(n3(mr))+", ms="+str(n3(ms)))

    if not verify_vm(s) or not verify_vm(r) or not verify_vm(I):
        raise NameError("Chosen V's and m's are not good")

    #------------------------- CapRIFping with Z1,Zr,Zs  -------------------------   
    if not (mq_ineq(s) and mq_ineq(r) and mq_ineq(I)):
        raise NameError("Impossible to cap since mq do not satisfy the needed inequalites")
    z1,zr,zs = zq(I),zq(r),zq(s)
    print("z1,zr,zs="+str((z1,zr,zs)))
    
    #  ------------------------- Edge potentials  -------------------------
    LQ={ (x,y): lq(x,y) for (x,y) in cartesian_product([[I,r,s],[I,r,s]]) }

    #  ------------------------- Find epsilon > 0  ------------------------- 
    epsmax = min(s/2, min([l for (l,q) in list(LQ.values())]))
    eps = find_epsilon(m1,mr,ms,QQ(n(epsmax)))
    print("Epsilon = "+str(n(eps)))

    all_curves()
    test_all()
    writing()
        
    return 0
    

