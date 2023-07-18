the_path = ""

# Main functions
load(the_path+"ternary_main_functions.sage")

# Conversions and simplifying
from sage.symbolic.expression_conversions import RingConverter
SR2RIF = RingConverter(RIF)

NR = lambda x: (RIF(x-1/2^54, x)).simplest_rational() 

#============================== Choosing values V,m,eps =============================

def choose_point(P):
    ''' A point between the vertex minimizing m1+mr+ms and the center'''
    weight = QQ(.9) # the weight of the found point
    v_minm1 = min(P.vertices(), key = lambda v: v[-3])
    v_minmr = min(P.vertices(), key = lambda v: v[-2])
    v_minms = min(P.vertices(), key = lambda v: v[-1])
    v_minzs =  max(P.vertices(), key = lambda v: v[-5])
    v_minm = [(x1+xr+4*xs)/6 for (x1,xr,xs) in zip(v_minm1, v_minmr, v_minms)] 
    center = P.center()
    the_point = [weight*v + (1-weight)*c for (v,c) in zip(v_minm, center)]
    return the_point

def verify_vm(q):
    '''since the polyhedron Pe1rs is an approximation, we still need to
    make sure that the choosen m1,mr,ms and vijk satisfy ALL the
    m-inequalities. '''
    print("    Checking m"+r_name(q)+"...")
    mq = mq_value(q)
    timer = time.time()
    for triangles in corona_generator(q):
        sum_v = sum(V(*t) for t in triangles)
        sum_angles = sum(tight_angle(t) for t in triangles)
        diff = RIF(abs(2*RIFpi - sum_angles))
        if not diff.contains_zero() and not mq >= - RIF(sum_v/diff):
            print(" !!! M"+str(r_name(q))+" = " + str(n(mq)) +  " < -sumv/diff = "+ str(n(- RIF(sum_v/diff)))+" for "+str([triangle_name(t) for t in triangles]))
            return False
    print("     Checking m"+r_name(q)+ " took " + str(time.time()-timer))
    return True


#============================== Find epsilon > 0 ==============================

def find_epsilon_rec(a, b, lb1, lbr, lbs):
    ''' Recurrence : true on a and false on b'''
    if b-a < QQ(0.0000000001):
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


def get_polyhedron(q, case):
    Pq = load(the_path+"output/"+str(case)+"_P_"+q+".sobj") #tmp
    #Pq = load(the_path+"output/"+str(case)+"_P_"+q+".sobj")
    print("  Found P"+q+": "+str(Pq))

    return Pq

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

zq = lambda q : QQ( max([  (-2 * RIFpi * V(*t)/tight_angle(t)).upper() for t in tripl(q) ])+1/2^53)


#============================== Testing on all triangles ==============================

def all_triangles():
    '''Returns a list of alll triangles we need to consider'''
    # an FM-triangle should at least contain half an s-disc -> lower bound on 16 times its squared area
    Smin=(16*((RIFpi/2)*s^2)^2).lower()

    T = []
    for (x,y,z) in ttriangles:
        if Q_inequality(x,y,z): # With it, we are allowed to consider only those with 1 contact
            T += [(y+z,RIF(x+z,x+z+2*s),RIF(x+y,x+y+2*s),x,y,z,False)]
            T += [(RIF(y+z,y+z+2*s),x+z,RIF(x+y,x+y+2*s),x,y,z,False)]*(r_name(x)!=r_name(y))
            T += [(RIF(y+z,y+z+2*s),RIF(x+z,x+z+2*s),x+y,x,y,z,False)]*(r_name(x)!=r_name(z) and r_name(y)!=r_name(z))
            T += [(0,RIF(x+z,x+z+2*s),RIF(x+y,x+y+2*s),x,y,z,True)]
            T += [(0, RIF(y+z,y+z+2*s), RIF(x+y,x+y+2*s),y,x,z,True)]*(r_name(x)!=r_name(y))
            T += [(0,RIF(x+z,x+z+2*s),RIF(y+z,y+z+2*s),z,y,x,True)]*(r_name(x)!=r_name(z) and r_name(y)!=r_name(z))
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

def test_triangle(a,b,c,ra,rb,rc,is_stretched):
    global nfeasable_num, epstight_num, good_num
    t_name = sorted_triangle_name((ra,rb,rc))
    if is_stretched: # r_a touches a: a = ()
        a = sqrt(b^2-ra^2)+sqrt(c^2-ra^2)
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
    if S >= 0 and not ((A*C).contains_zero() and (B.contains_zero() or (4*A*C/B^2).upper()>.78)):
        R=radiusABCD(a,b,c,ra,rb,rc,A,B,C,D)
        if R > s:    # not feasible because not saturated
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
            raise NameError("(case "+str(CASE)+") "+t_name+("_stretched" if is_stretched else " ")+" E<U: "+str(n(ee))+" < "+str(uu)+" with "+infos)

    # The actual precision does not allow to conclude -> refine
    # taking into account that some values out of a,b,c are not really
    # intervals: we should not divide them 
    EPS = 1e-9 # Took it randomly...
    
    glue_ends = lambda point: [point] if point.diameter()<=EPS else [RIF(z) for z in list(RIF(point).bisection())]
    if is_stretched:
        a = RIF(0)
    new_triples = cartesian_product([glue_ends(x) for x in [a,b,c]])
    
    if len(new_triples) == 1: #all of them are points: we are blocked here
        save_triangle(a,b,c,ra,rb,rc,"_blocked")
        raise NameError("(case "+str(CASE)+") "+"We're blocked on: " +t_name+("_stretched" if is_stretched else " ")+" E<U: "+" "+str((a,b,c))) 

    for (a2,b2,c2) in new_triples:
        test_triangle(a2,b2,c2,ra,rb,rc, is_stretched)
    return 0
    
    
def test_all():
    T = all_triangles()
    for (a,b,c,ra,rb,rc,is_stretched) in T:
        print(triangle_name((ra,rb,rc))+(" stretched " if is_stretched else " ")+": testing...")
        test_triangle(a,b,c,ra,rb,rc,is_stretched)
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

def testing_cases_from_list(l):
    for c in l:
        load('ternary_general.sage')
        with open("testing_results.txt", "a") as myfile:
            myfile.write(str(c))
            try:
                main(c)
                all_curves()
                test_all()
                myfile.write(" OK\n")
            except NameError as err:
                print(err)
                print("Got an error proceeding this case :-(")
                myfile.write(" not OK"+str(err)+"\n")
        reset()
    print("DONE")

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
   
    
def do_things(case):
    global V2, m1, mr, ms, eps, z1, zr, zs, LQ
    init(case)
    
    # ------------------------- Computing Vxyz  -------------------------
    M,v=matrix_vector(coronas)
    V2=M.inverse()*v

    #  ------------------------- Polyhedron  -------------------------
    Pvmeps = get_polyhedron("e1rs",case)

    #  ------------------------- Choose a point Vm -------------------------
    vm_point = choose_point(Pvmeps) # 7 free variables

    m1,mr,ms = map(NR, vm_point[-3:])
    
    v_dict = {v:val for v,val in zip([vr1r,vs1s,v1r1,vsrs,v1s1,vrsr], vm_point[:-3])}
    V2 = V2.subs(v_dict) # get rid of variables
    V2 = [NR(QQ(RIF(x).center())) for x in V2]
    print(list(map((lambda x: (x,n3(v_dict[x]))),v_dict)))
     
    print(" m1="+str(n3(m1))+", mr="+str(n3(mr))+", ms="+str(n3(ms)))

    #if not verify_vm(s) or not verify_vm(r) or not verify_vm(I):
    #    raise NameError("Chosen V's and m's are not good")

    #------------------------- CapRIFping with Z1,Zr,Zs  -------------------------   
    if not (mq_ineq(s) and mq_ineq(r) and mq_ineq(I)):
        raise NameError("Impossible to cap since mq do not satisfy the needed inequalites")
    z1,zr,zs = zq(I),zq(r),zq(s)
    print("z1,zr,zs="+str((z1,zr,zs)))
    
    #  ------------------------- Edge potentials  -------------------------
    LQ={ (x,y): lq(x,y) for (x,y) in cartesian_product([[I,r,s],[I,r,s]]) }

    #  ------------------------- Find epsilon > 0  ------------------------- 
    epsmax = min(s/2, min([l for (l,q) in list(LQ.values())]))
    eps = NR(find_epsilon(m1,mr,ms,QQ(n(epsmax))))
    if not test_eps(m1,mr,ms,eps):
        raise NameError("Epsilon rounding NR was not precise enough")
    
    print("Epsilon = "+str(n(eps)))


    #save_curves_for_paper()
    #writing()
    test_all()

        
    return 0
    

#proved_cases_a = [53,54,55,56,66,76,77,79,93,108,115,116,118,129,131, 146]
#TODO   131,146
