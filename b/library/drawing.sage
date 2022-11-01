# independent of dics number

#============================== Draw Curves ============================================================

#------------------------------ Thomas' way to draw ------------------------------

# to draw lines instead of plot (which is capricious)
#seg = lambda a,b,k: [RDF(a+i*(b-a)/k) for i in range(k+1)]

seg = lambda a,b,k: [QQ(n(a+i*(b-a)/k)) for i in range(k+1)]

# for each edge xy, study triangles xyz with contact xz and yz and length of edge xy varying from tight to stretched


var('x')

def curves_for_paper(ra,rb,rc):
    ''' Edge b (between ra-center and rb-center )varies from contact to stretched '''
    bmin = ra+rc
    bmax = min(s/10+sqrt((rb+rc)**2-rb**2) + sqrt((rb+ra)**2-rb**2), ra+rb+rb+rc-QQ(.0000001))
    Ep = line([[x, E(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],thickness = .5)
    Up = line([[x, Unoz(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],color='red', linestyle='dotted', thickness = .5)
    U2p = line([[x, U2(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,200)],color='red', thickness = 1.3)
    tD = sign_change_point(lambda x: dist(x,rb+rc,ra+rb,rb,ra,rc,radius(rb+rc,x,ra+rb,ra,rb,rc)), ra+rc, bmax, QQ(.05))
    Dp=line([[x, dist(x,rb+rc,ra+rb,rb,ra,rc,radius(rb+rc,x,ra+rb,ra,rb,rc))] for x in seg(tD-.05,tD+.05,50)],color='magenta',thickness = .5)
    Zp = line([[x, U(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],color = 'red',linestyle='dashed', thickness = .8)
   # Zp = plot(z(rb), color = 'red', xmin=ra+rc, xmax=bmax, thickness = .8)
    #Epsp = line([(ra+rc+eps,0), (ra+rc+eps,.1)], linestyle='dashed', color = "orange", thickness = .5)
    L =  line([(LQ[(ra,rc)][0],-0.03), (LQ[(ra,rc)][0],0.03)], color = "green", thickness = 1)
    return U2p+Up+Ep+Dp+L+Zp


def curves(ra,rb,rc):
    ''' Edge b (between ra-center and rb-center )varies from contact to stretched '''
    #print(ra,rb,rc)
    bmin = ra+rc
    bmax = min(s/10+sqrt((rb+rc)**2-rb**2) + sqrt((rb+ra)**2-rb**2), ra+rb+rb+rc-QQ(.0000001))
    Ep = line([[x, E(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],thickness = .5)
    Up = line([[x, Unoz(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],color='red',linestyle='dotted', thickness = .5)
    U2p = line([[x, U2(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,200)],color='brown',linestyle='dashed', thickness = 1)
    tD = sign_change_point(lambda x: dist(x,rb+rc,ra+rb,rb,ra,rc,radius(rb+rc,x,ra+rb,ra,rb,rc)), ra+rc, bmax, QQ(.05))
    Dp=line([[x, dist(x,rb+rc,ra+rb,rb,ra,rc,radius(rb+rc,x,ra+rb,ra,rb,rc))] for x in seg(tD-QQ(.05),tD+QQ(.05),50)],color='magenta',thickness = .5)
    Zp = line([[x, U(rb+rc,x,ra+rb,ra,rb,rc)] for x in seg(bmin,bmax,100)],color = 'red', thickness = .8)
    #Zp = plot(z(rb), color = 'red', xmin=ra+rc, xmax=bmax, thickness = .8)
    Epsp = line([(ra+rc+eps,QQ(0)), (ra+rc+eps,QQ(.1))], linestyle='dashed', color = "orange", thickness = .5)
    L = line([(LQ[(ra,rc)][0],QQ(-0.03)), (LQ[(ra,rc)][0],QQ(0.03))], color = "green", thickness = 1)
    return Zp+U2p+Up+Ep+Dp+Epsp+L

def all_curves(name=""):
    rttriangles = tripl(I)+tripl(r)+tripl(s)
    for (ra,rb,rc) in rttriangles:
        (curves(ra,rb,rc)).save(dpi=500, filename = "curves/"+str(CASE)+"_"+r_name(ra)+r_name(rb)+r_name(rc)+name+".pdf")



def save_curves_for_paper(name=""):
    rttriangles = tripl(I)+tripl(r)+tripl(s)
    for (ra,rb,rc) in rttriangles:
        (curves_for_paper(ra,rb,rc)).save(dpi=500, filename = "curves/paper_"+str(CASE)+"_"+r_name(ra)+rname(rb)+rname(rc)+name+".pdf") 

#======================================== Draw triangles ========================================

def save_triangle(a,b,c,ra,rb,rc, filename):
    '''Gets a figure and saves its as filename.pdf
    '''
    f = draw_triangle(a,b,c,ra,rb,rc)
    aa,bb,cc = numerical_approx(a,digits=5),numerical_approx(b,digits=5),numerical_approx(c,digits=5)
    (f.plot()).save('figures/'+str(CASE)+"_"+r_name(ra)+r_name(rb)+r_name(rc)+"_"+str(aa)+"-"+str(bb)+"-"+str(cc)+"_"+filename+'.pdf',dpi=100)

dcol = lambda x : (.9,.8,0) if x>r else (.4,.55,.9) if x>s else (.6,0,0)
    
def draw_triangle(a,b,c,ra,rb,rc):
    '''Returns a figure representing the triangle+its discs
    '''
    figure = Graphics()
    cosc = (a**2+b**2-c**2)/(2*a*b)
    sinc = sqrt(1-cosc**2)
    A,B,C = (b*cosc, b*sinc), (a,0), (0,0),
    figure += circle(A, ra, fill=False, edgecolor = dcol(ra),zorder=0)
    figure += circle(B, rb, fill=False, edgecolor = dcol(rb),zorder=0)
    figure += circle(C, rc, fill=False, edgecolor = dcol(rc),zorder=0)
    figure += polygon([A,B,C], rgbcolor='black', fill=False, zorder=2)
    return figure
    
def draw_triangle_with_info(a,b,c,ra,rb,rc):     
    # V(ra,RC,rb)
    figure = draw_triangle(a,b,c,ra,rb,rc)
    figure+=text("V"+r_name(ra)+r_name(rc)+r_name(rb)+" = "+str(numerical_approx(V(ra,rc,rb),digits=3)), (-1.1,-1.1), rgbcolor='gray',zorder=2)
    # Uv(ra,RC,rb)
    figure+=text("Uv = "+str(numerical_approx(Uv(a,c,b,ra,rc,rb),digits=3)), (-1.1,-1.3), rgbcolor='black',zorder=2)
    # V(ra,RB,rc)
    figure+=text("V"+r_name(ra)+r_name(rb)+r_name(rc)+" = "+str(numerical_approx(V(ra,rb,rc),digits=3)), (a+1.1,-1.1), rgbcolor='gray',zorder=2)
    # Uv(ra,RB,rc)
    figure+=text("Uv = "+str(numerical_approx(Uv(a,b,c,ra,rb,rc),digits=3)), (a+1.1,-1.3), rgbcolor='black',zorder=2)
    # V(rb,RA,rc)
    figure+=text("V"+r_name(rb)+r_name(ra)+r_name(rc)+" = "+str(numerical_approx(V(rb,ra,rc),digits=3)), (b*cosc+1.1, b*sinc+1.3), rgbcolor='gray',zorder=2)
    # Uv(rb,RA,rc)
    figure+=text("Uv = "+str(numerical_approx(Uv(b,a,c,rb,ra,rc),digits=3)), (b*cosc+1.1, b*sinc+1.1), rgbcolor='black',zorder=2)

    # U(triangle), E(triangle)
    figure+=text("Uv(T) = "+str(numerical_approx(U(a,b,c,ra,rb,rc),digits=5)), (a/2,-1.3), rgbcolor='red',zorder=2)
    figure+=text("U2(T) = "+str(numerical_approx(U2(a,b,c,ra,rb,rc),digits=5)), (a/2,-1.5), rgbcolor='red',zorder=2)
    figure+=text("E(T) = "+str(numerical_approx(E(a,b,c,ra,rb,rc),digits=5)), (a/2,-1.7), rgbcolor='blue',zorder=2)

    figure.axes(False)
    return figure


#======================================== Draw Corona ========================================

# turn (x,y) by 'a' counterclockwise
turn = lambda x,y,a: (x*cos(a)-y*sin(a), x*sin(a)+y*cos(a))

def find_third_center(r1,x1,y1, r2,x2,y2, r3):
    '''Returns the coordinates of a new disc with two neighbors given  '''
    dx, dy = x2-x1, y2-y1
    a1, ratio = angle(r3+r2,r3+r1,r1+r2), (r1+r3)/(r1+r2)
    xv, yv = turn(ratio*dx, ratio*dy, a1)
    center = (x1+xv, y1+yv)
    return(center)

def draw_corona(r_0, corona_s):
    corona = list(map(choose_r, corona_s))
    figure = Graphics()
    figure += circle((0,0), r_0, fill=True, edgecolor = 'black', facecolor = dcol(r_0),zorder=0)
    x,y,rr = r_0+corona[0], 0, corona[0]
    x1,y1,r1 = x,y,rr
    figure += circle((x,y), rr, fill=True, edgecolor = 'black', facecolor = dcol(rr),zorder=0)
    for t in corona[1:]:
        xn,yn = find_third_center(r_0,0,0, rr,x,y, t)
        figure += circle((xn,yn), t, fill=True, edgecolor = 'black', facecolor = dcol(t),zorder=0)
        # V(ra,RC,rb)
        #figure+=text("V"+r_name(rr)+r_name(r_0)+r_name(t)+" = "+str(numerical_approx(V(rr,r_0,t),digits=3)), (1.1*(x+xn),1.1*(y+yn)), rgbcolor='gray',zorder=2)
        figure += polygon([(x,y), (xn,yn), (0,0)], rgbcolor = 'black', fill=False, zorder=2 )
        rr,x,y=t,xn,yn
    # the last with the first
    #figure+=text("V"+r_name(rr)+r_name(r_0)+r_name(r1)+" = "+str(numerical_approx(V(rr,r_0,r1),digits=3)), (1.1*(x+x1),1.1*(y+y1)), rgbcolor='gray',zorder=2)
    figure += polygon([(x,y), (x1,y1), (0,0)], rgbcolor = 'black', fill=False, zorder=2 )
    # figure.axes(False)
    return figure



# ------------------------------------------------ Drawing packings ------------------------


def find_center(r, neighbors):
    '''Returns the coordinates of a new disc with two neighbors given as a
    pair of dictionaries  '''
    (x1,y1), (x2,y2) = neighbors[0]['xy'], neighbors[1]['xy']
    dx, dy = x2-x1, y2-y1
    r1, r2  =  neighbors[0]['r'], neighbors[1]['r']
    a1, ratio = angle(r+r2,r+r1,r1+r2), (r1+r)/(r1+r2)
    xv, yv = turn(ratio*dx, ratio*dy, a1)
    center = (x1+xv, y1+yv)
    return(center)

def draw_packing(d,discs,shifts,a_range):
    '''Gets a list of discts describing the packing and returns a figure
    representing the domain of this packing'''
    
    xmin,xmax,ymin,ymax = tuple(a_range)
    col = lambda x: (227/256,15/256,15/256) if x=='s' else (255/256,212/256,40/256) if x=='1' else (84/256,151/256,255/256)
    figure = Graphics()
    for shift in shifts:
        for disc in discs:
            figure += circle(vector(disc['xy'])+vector(shift), disc['r'], fill=True, facecolor = col(disc['rname']), edgecolor = 'black', thickness = .1, zorder=1)

    R,r = max([disc['r'] for disc in discs]) , min([disc['r'] for disc in discs])
    Rname,rname = min([disc['rname'] for disc in discs]), max([disc['rname'] for disc in discs])
    rr = [x for x in [disc['r'] for disc in discs if disc['rname']=='r']][0]
    border_width = (xmax-xmin)/10
    white_rec = polygon([(xmin-border_width,ymin),(xmax+border_width,ymin),(xmax+border_width,ymin-border_width*2),(xmin-border_width,ymin-border_width*2)],color='white',zorder=1)
    figure+=white_rec
    figure += text("d="+str(n(d,digits = 8)),((xmin+xmax)/2,ymin-border_width/3),horizontal_alignment = 'right',color='black',fontsize='large' )
    figure += text(rname+"="+str(n(r,digits=8)), ((3*xmin+4*xmax)/7,ymin-border_width/3),horizontal_alignment = 'left', color=col(rname),fontsize='large' )
    figure += text("r="+str(n(rr,digits=8)), ((3*xmin+4*xmax)/7,ymin-6*border_width/7), horizontal_alignment = 'left', color=col('r'),fontsize='large' )
    if Rname!='1':
        figure += text(Rname+"="+str(n(R,digits=8)), ((3*xmin+4*xmax)/7,ymin-6*border_width/7), horizontal_alignment = 'left', color=col(Rname),fontsize='large' )
    figure.set_axes_range(xmin,xmax,ymin-border_width,ymax)
    return figure


#============================== Draw Triangulated Ternary Packings ============================================================


def encoding(case,r,s):
    '''Returns the encoding of the case case where the radii are equal to
    r and s'''
    from packings3discs import packings
    return (packings(r,s))[case-1]

def halfplane_sign(p, a, b):
    '''Returns the sign of the point p with respect to the line passing
    by a and b'''
    px,py = p
    ax,ay = a
    bx,by = b
    # y = mx + n  is the line of the points a and b
    m = (by-ay)/(bx-ax)
    n = ay-m*ax
    return (py - (m*px + n)).sign()


def find_center(r, neighbors):
    '''Returns the coordinates of a new disc with given neighbors as a
    pair of dictionaries'''
    (x1,y1), (x2,y2) = neighbors[0]['xy'], neighbors[1]['xy']
    dx, dy = x2-x1, y2-y1
    r1, r2  =  neighbors[0]['r'], neighbors[1]['r']
    a1, ratio = angle(r+r2,r+r1,r1+r2), (r1+r)/(r1+r2)
    xv, yv = turn(ratio*dx, ratio*dy, a1)
    center = (x1+xv, y1+yv)
    return(center)

def touch(r1, x1, y1, r2, x2, y2, epsilon):
    return sqrt((x1-x2)^2 + (y1-y2)^2) - (r1+r2) < epsilon

# TODO: test it better
def other_neighbors(r, xy, neighbors, discs, epsilon):
    '''Returns the neighbors of a new disc except the already given ones
    with given precision epsilon '''
    on = []
    n_discs = [d for d in discs if not d in neighbors]
    for disc in n_discs :
        if touch(r, xy[0], xy[1], disc['r'], disc['xy'][0], disc['xy'][1], epsilon):
            on.append(discs.index(disc))
    #if len(on) > 0:
    #    print("Other neighbors: "+str(len(on)))
    return(on)

# TODO: test better ?
def find_e(r,s):
    '''Returns the mimimal distanc between non-tangent discs in the
    compact packing'''
    radii = [1,r,s]
    pairs = [(r1,r2) for r1 in radii for r2 in radii if r1<=r2]
    min_dist = 2
    for (r_top,r_bottom) in pairs:
        for (r_left, r_right) in pairs:
            a_lt = angle(r_left+r_bottom, r_left+r_top, r_top+r_bottom)
            a_rt = angle(r_left+r_bottom, r_left+r_top, r_top+r_bottom)
            dist = sin(a_lt)*(r_left+r_top) + sin(a_rt)*(r_right+r_top) - r_left -r_right
            if dist < min_dist:
                min_dist = dist
    return min_dist

# turn (x,y) by 'a' counterclockwise
turn = lambda x,y,a: (x*cos(a)-y*sin(a), x*sin(a)+y*cos(a))

# OK
def discs(case,r,s):    
    '''Returns the geometric properties of the packing'''
    l = (encoding(case,r,s))[1]
    name = lambda x : '1' if x>r else 'r' if x>s else 's'
    disc_0 = {'r':RIF(l[0]), 'xy':(0,0),'rname':name(l[0])}
    disc_1 = {'r':RIF(l[1]), 'xy':(disc_0['r']+RIF(l[1]),0), 'rname':name(l[1])}
    discs = [disc_0,disc_1]
    # the precision: two discs can not be at distance < epsilon and do not touch
    epsilon = find_e(r,s)/2
    for entry in l[2:]:
        R = RIF(entry[2])
        neighbors = list(entry[:2])
        xy = find_center(R, [discs[i] for i in neighbors])
        d = {'r':R, 'xy':xy, 'rname':name(R)}
        index = l.index(entry)
        if len(entry)==4:
            twin = entry[3]
            (discs[twin]).update({'twin':index})
            d.update({'twin':twin})
        discs.append(d)
        #print("Disc number "+str(index))# + ":   " + str([r,xy,neighbors[:2]]))


    periods = []
    for disc in discs:
        if 'twin' in disc.keys():
            disctwin = discs[disc['twin']]
            x1,y1 = disc['xy']
            x2,y2 = disctwin['xy']
            v = vector((x1-x2,y1-y2))
            if (not v in periods) and (not (-v) in periods) :
                periods.append(v)
    return discs,periods



def save_ter(case):
    all_discs,periods = discs(case,r,s)
    shifts = [i*v1+j*v2 for i in [-2,-1,0,1,2] for j in [-2,-1,0,1,2] for v1 in periods for v2 in periods]
    f = draw_packing(d_opt,all_discs,shifts,[-4,4,-4,4])
    f.save(dpi=100, axes = False, filename = "packings_drawings/"+str(case)+".pdf")
