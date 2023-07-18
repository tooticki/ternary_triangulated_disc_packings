the_path=''

# Main functions
load(the_path+"ternary_main_functions.sage")

#==============================  Global Variables ==============================

from sage.symbolic.expression_conversions import RingConverter
SR2RIF = RingConverter(RIF)

ptype = 'e1rs'

#============================== Epsilon polyhedron  ==============================

def eps_polyhedron(eps):
    ''' Returns the polyhedron on V(all R6), m1,mr,ms st dE > dU with eps, not precise !!! (RDF) '''
    eqs = [(0,1,0,0),(0,0,1,0),(0,0,0,1)]  # m1,mr,ms >= 0
    var('a','b','c', 'M1', 'Mr', 'Ms')
    endps = lambda  x: [0] if RIF(x)==0 else RIF(x).endpoints()
    for (ra,rb,rc) in ttriangles:
        for x in [a,b,c]:
            De = dE(ra,rb,rc,eps,x)
            le = De.endpoints()
            Du = dU(ra,rb,rc,eps, x)
            lu1 = endps(Du.coefficient(M1))
            lur = endps(Du.coefficient(Mr))
            lus = endps(Du.coefficient(Ms))
            endpoints = [(de,du1,dur,dus) for de in le for du1 in lu1 for dur in lur for dus in lus]
            # The intervals of de,du are not precise, so we check all combinations of endpoints
            for (de,du1,dur,dus) in endpoints:
                eqs += [(QQ(de), -QQ(du1), -QQ(dur), -QQ(dus))]
    eqs9D = [ [e[0]] + [0,0,0,0,0,0] + list(e[1:]) for e in eqs]
    return Polyhedron(ieqs = eqs9D)

#============================== Vm polyhedron  ==============================

def Vm_polyhedron(q):
    '''Given q, returns a Q^9- polyhedron[ vr1r,vsrs,etc]+[m_1, m_r, m_s]
    taking into account the lower bound on m_q ensuring potential
    positivity around any vertex in a packing (consider only
    "important" configs around vertices: should be checked with all
    config afterwards by the function verify_m(q,mq))

    '''
    timer = time.time()
    m_q_vector = [1,0,0] if q>r else [0,1,0] if q>s else [0,0,1]    
    qq = lambda x : RIF(x).simplest_rational() # simplest rational...
    
    exact_coronas_num = 0
    Rifr.<vr1r, vs1s, v1r1, vsrs, v1s1, vrsr> = PolynomialRing(RIF)
    vlist =  [vr1r, vs1s, v1r1, vsrs, v1s1, vrsr]
    ineqs = []
    l = list(corona_generator(q)) 
    print("Polyhedron Vm" + r_name(q) + ": we have " + str(len(l)) + " coronas to check")
    counter,t,t0=0, time.time(), time.time()
    for triangles in l:
        counter+=1
        if time.time()-t > 60:
            t = time.time()
            print("  " + str(t-t0) + " seconds of computation: we have treated " + str(counter) + " coronas")
        sum_v = sum(V(*t) for t in triangles)
        sum_angles = sum(tight_angle(t) for t in triangles)
        diff = abs(2*RIFpi - sum_angles)

        if diff.contains_zero():
            exact_coronas_num += 1
        else:
            lb = -sum_v / diff
            lbp = Rifr(lb)
            # mq >= lbp -> -lbp+mq >= 0
            ineq = list(map(qq, [-lbp.constant_coefficient()]+[-lbp.coefficient({y:1}) for y in vlist]+m_q_vector))
            ineqs+=[ineq]

    print("Eqs done, making the polyhedron")
    PP = Polyhedron(ieqs = ineqs)

    print("     Polyhedron for m"+str(r_name(q))+ " took " + str(time.time()-timer)+" | "+str(PP))
    if exact_coronas_num > 3:
            print ("     Too many exact coronas around " + r_name(q))  
    return PP


#============================== Treating types of polyhedra  ==============================

def get_all_polyhedra(case):
    Pe1rs = get_one_polyhedron("e1rs",case)
    if Pe1rs.is_empty():
        raise NameError("Pe1rs is empty")
    if not Pe1rs.is_compact():
        raise NameError(str(case)+" Pe1rs is not compact")
    return Pe1rs
        

def get_one_polyhedron(q, case):
    ''' If a polyhedron of type string(q) was already saved for this case, load it
    else, compute it and save it'''
    Pq = 0
    try:        
        Pq = load(the_path+"output/"+str(case)+"_P_"+q+".sobj") #TODO upper
        print("  Found P"+q+": "+str(Pq))
    except FileNotFoundError:
        print("  Making P"+q+": "+str(Pq))
        if q=="e":
            Pq = eps_polyhedron(0)
            if Pq.is_empty():
                raise NameError("Pe is empty for eps=0")
        elif q=="1":
            Pq = Vm_polyhedron(I)
        elif q=="r":
            Pq = Vm_polyhedron(r)
        elif q=="s":
            Pq = Vm_polyhedron(s)
        elif q=="e1":
            Pq = get_one_polyhedron("e",case).intersection(get_one_polyhedron("1",case))
        elif q=="rs":
            Pq = get_one_polyhedron("r",case).intersection(get_one_polyhedron("s",case))
        elif q=="e1rs":
            Pq = get_one_polyhedron("e1",case).intersection(get_one_polyhedron("rs",case))
        
        print(str(case)+"  P"+q+" was computed: "+str(Pq))
        save(Pq, the_path+"output/"+str(case)+"_P_"+q+".sobj")#TODO add output
        print(str(case)+"  P"+q+" was saved")
    return Pq
    

        
def main(case, ptype):
    global V2, m1, mr, ms, eps, z1, zr, zs, LQ
    init(case)
    
    # Vxyz
    M,v=matrix_vector(coronas)
    V2=M.inverse()*v

    # Polyhedra
    Pvmeps = get_one_polyhedron(ptype, case)
        
    return 0
    
