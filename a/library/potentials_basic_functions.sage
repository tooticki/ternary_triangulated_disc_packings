# independent of dics number

#------------------------------ Finally, vertex potentials ------------------------------


def Uv(a,b,c,ra,rb,rc):
    '''Uv where v is opposite to the edge b (center of a circle rb) in the
    triangle on edges a,b,c st the circle of radius rx is opposite to x '''
    return min(zq_value(rb), V(ra,rb,rc) + mq_value(rb) * abs(angle(b,c,a) - tight_angle((ra,rb,rc))))

def Uvnoz(a,b,c,ra,rb,rc):
    '''Uv where v is opposite to the edge b (center of a circle rb) in the
    triangle on edges a,b,c st the circle of radius rx is opposite to x '''
    return (V(ra,rb,rc) + mq_value(rb) * abs(angle(b,c,a) - tight_angle((ra,rb,rc))))


def U(a,b,c,ra,rb,rc):
    ''' The sum of 3 vertex potentials '''
    return (Uv(a,b,c,ra,rb,rc) + Uv(b,c,a,rb,rc,ra) + Uv(c,a,b,rc,ra,rb))

def Unoz(a,b,c,ra,rb,rc):
    ''' The sum of 3 vertex potentials '''
    return (Uvnoz(a,b,c,ra,rb,rc) + Uvnoz(b,c,a,rb,rc,ra) + Uvnoz(c,a,b,rc,ra,rb))




#------------------------------ Potentials and excess derivatives ------------------------------

def dE(Ra,Rb,Rc,eps,x):
    '''Returns the interval RIF(dE/dx) where x is on of the sides of a
    triangle: a,b or c '''
    var('a','b','c', 'ra', 'rb', 'rc')
    # to substitute edge names with edge length intervals
    edges={a:(RIF(Rb+Rc,Rb+Rc+eps)),b:(RIF(Ra+Rc,Ra+Rc+eps)),c:(RIF(Ra+Rb,Ra+Rb+eps))}
    de = derivative(symb_E(a,b,c,Ra,Rb,Rc), x)
    return RIF(de.subs(edges))


def dU(Ra,Rb,Rc,eps,x):
    ''' Returns the interval RIF(dU/dx) where x is on of the sides of a
    triangle: a, b or c '''
    var('a','b','c', 'ra', 'rb', 'rc', 'M1', 'Mr', 'Ms')
    mq_var = lambda x: M1 if x>r else Mr if x>s else Ms
    # to substitute edge names with edge length intervalsR
    edges={a:(RIF(Rb+Rc,Rb+Rc+eps)),b:(RIF(Ra+Rc,Ra+Rc+eps)),c:(RIF(Ra+Rb,Ra+Rb+eps))}
    # sum over the three vertices
    ca=RIF(abs(derivative(symb_angle(a,b,c),x).subs(edges)))
    cb=RIF(abs(derivative(symb_angle(b,c,a),x).subs(edges)))
    cc=RIF(abs(derivative(symb_angle(c,a,b),x).subs(edges)))
    if is_infinity(ca) or is_infinity(cb) or is_infinity(cc):
        raise NameError("One of dU coefficients is = infty : " + str((Ra,Rb,Rc,eps, x)) + " ca,cb,cc: "+str((ca,cb,cc)))
    z = ca*mq_var(Ra) + cb*mq_var(Rb) +  cc*mq_var(Rc)
    return SR(z)

##########################################
# add edge potential to vertex potential #
##########################################

def U2(a,b,c,ra,rb,rc):
    R=radius(a,b,c,ra,rb,rc)
    #if RIF(R).diameter > .01:
    #    print (list(map(n, RIF(R).endpoints())))
    u = U2withR(a,b,c,ra,rb,rc,R)
    return u

def U2withR(a,b,c,ra,rb,rc,R):
    #if RIF(R).diameter > .01:
    #    print (list(map(n, RIF(R).endpoints())))
    u = U(a,b,c,ra,rb,rc)
    
    (la,qa) = LQ[(rb,rc)]
    if a > la:
        u +=qa*dist(a,b,c,ra,rb,rc,R)
        
    (lb,qb) = LQ[(ra,rc)]
    if b > lb:
        u += qb*dist(b,c,a,rb,rc,ra,R)

    (lc,qc) = LQ[(ra,rb)]
    if c > lc:
        u += qc*dist(c,a,b,rc,ra,rb,R)
    return u
