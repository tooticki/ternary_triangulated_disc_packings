thepath = ''
# Input from the table
attach(thepath+'library/input.sage')

# All calculatory functions are here
attach(thepath+'library/small_functions.sage')

# All drawing functions are here
attach(thepath+'library/drawing.sage')


Qpi = QQ(n(pi))
#########################################################################################
# Binary compact packings
#########################################################################################

var('r')

bin_r=[
AA(ZZ[r](r^4-10*r^2-8*r+9).roots(AA)[0][0]),
AA(ZZ[r](r^8 - 8*r^7 - 44*r^6 - 232*r^5 - 482*r^4 - 24*r^3 + 388*r^2 - 120*r + 9).roots(AA)[2][0]),
AA(ZZ[r](8*r^3+3*r^2-2*r-1).roots(AA)[0][0]),
AA(sqrt(2)-1), 
AA((2*sqrt(3)+1-2*sqrt(1+sqrt(3)))/3),
AA(sin(pi/12)/(1-sin(pi/12))),
AA((sqrt(17)-3)/4),
AA(2/sqrt(3)-1),
AA(5-2*sqrt(6))]


# corresponding density/pi (algebraic value)
bin_d=[
ZZ[r](27*r^4 + 112*r^3 + 62*r^2 + 72*r - 29).roots(AA)[1][0],
ZZ[r](r^8 - 4590*r^6 - 82440*r^5 + 486999*r^4 - 1938708*r^3 + 2158839*r^2 - 1312200*r + 243081).roots(AA)[2][0],
ZZ[r](1024*r^3 - 692*r^2 + 448*r - 97).roots(AA)[0][0],
ZZ[r](2*r^2 - 4*r + 1).roots(AA)[0][0],
ZZ[r](944784*r^4 - 3919104*r^3 - 2191320*r^2 - 1632960*r + 757681).roots(AA)[0][0],
ZZ[r](144*r^4 + 9216*r^3 + 133224*r^2 - 127104*r + 25633).roots(AA)[2][0],
ZZ[r](4096*r^4 + 2924*r^2 - 289).roots(AA)[1][0],
ZZ[r](108*r^2 + 288*r - 97).roots(AA)[1][0],
ZZ[r](144*r^4 - 4162200*r^2 + 390625).roots(AA)[2][0]]

bin_names = ["b"+str(i+1) for i in range(9)]

#Contains  packing_name : (radius, density)
bin_dict = {"b"+str(i+1):  (RIF(bin_r[i]),RIF(bin_d[i]*pi))  for i in range(9)}
for x in bin_names:
    print(x,bin_dict[x])
        
    
#########################################################################################
# Ternary compact packings
#########################################################################################

ter_names = [i+1 for i in range(164)]
ter_dict = {}

def ternary_init():
    global ter_dict
    for case in ter_names:
        (r,s, disc_nums_1rs, triangles_nums, coronas, dens) = input_values(case)
        ter_dict.update({case:(r,s,s/r,RIF(dens))})


#########################################################################################
# Candidates for counter-examples
#########################################################################################

eps = QQ(0.2)

def ternary_with_similar_ratio(bin_name):
    ratio, density   = bin_dict[bin_name]
    lr,ls,lsr = [],[],[]
    for case in ter_names:
        ter_r, ter_s, ter_sr, ter_d  = ter_dict[case]
        if ter_d < density:
            if abs(ter_r-ratio) < eps:
                lr+=[{'case':case, 'r':ter_r, 'R':1, 'ter_d':ter_d}]
            if abs(ter_s-ratio) < eps:
                ls+=[{'case':case, 'r':ter_s, 'R':1, 'ter_d':ter_d}]
            if abs(ter_sr-ratio) < eps:
                lsr+=[{'case':case, 'r':ter_s, 'R':ter_r, 'ter_d':ter_d}]
    return(('r1',lr),('s1',ls),('sr',lsr))
           
           
def denser_than_any_bin():
    l = []
    dmaxbin = max(bin_dict[bname][1] for bname in bin_names)
    for case in ter_names:
        if ter_dict[case][-1] > dmaxbin:
            l.append(case)
    return l

# r8 <-> r7
d0=(0,0,1)
d1=(2,0,1)
d2=stick(d0,d1,r)
d3=stick(d0,d2,1)
# r8: when d1 touches d3, i.e., d3[0]=1
r8=[x for x in find_roots(d3[0]-1,0,1) if x >0 and x<1][0]
denser_than_bin = [23, 26, 32, 35, 36, 38, 42]
ter_names = [i for i in ter_names if not i in denser_than_bin]

#########################################################################################
# Candidates density
#########################################################################################

def base_area(bin_name):
    R,r = 1,bin_dict[bin_name][0]
    if bin_name == "b1":
        return 6*d_area(R,R,r)+2*d_area(R,r,r)
    elif bin_name == "b2":
        return 6*d_area(R,R,r)+3*d_area(R,r,r) + d_area(r,r,r)  + 3*d_area(R,R,R)
    elif bin_name == "b4":
        return 4*R**2
    elif bin_name == "b7":
        return  4*d_area(R,R,r) + 2*d_area(R,r,r)
    elif bin_name == "b8":
        return area(2*R, 2*R, 2*R)
    elif bin_name == "b9":
        return area(2*R, 2*R, 2*R)
    else:
        return 0

def almost_compact_density_r_smaller(bin_name, r, R):
    tile_area = base_area(bin_name) * (R**2 if bin_name!="b6" else (r/(bin_dict["b6"][0]))**2)
    if bin_name == "b1":
        return 2*Qpi*(R**2+r**2) / tile_area
    elif bin_name == "b2":
        return Qpi*(2*R**2+3*r**2) / tile_area
    elif bin_name == "b3":
        return Qpi*(R**2+2*r**2) / (2*R * (R + r + sqrt((R+r)**2-R**2)))
    elif bin_name == "b4":
        a = tight_angle((R,R,r))
        if sin(a) < QQ(.5): # r is too small
            return 0
        return (Qpi*(R**2 + r**2) / ( 8*sin(a)*cos(a)*R**2 ))
    elif bin_name == "b5":
        return (2*Qpi*R**2+7*Qpi*r**2)/6/sqrt(3)/R**2
    elif bin_name == "b6":
        a = (2*Qpi - tight_angle((r,r,r)) - 2*tight_angle((R,r,r)))/2
        if RIF(R)  >  RIF((r+R)*sin(a)): # r is too small
            return 0
        return Qpi*(3*r**2 + QQ(.5)*R**2) / (3*d_area(R,r,r)+ d_area(r,r,r) + 3*sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b7":
        a = Qpi-tight_angle((R,r,R))
        if RIF((3*r)/(r+R)) < RIF(cos(a)): # r is too small
            return 0               
        return Qpi*(R**2 + 2*r**2) / (4*d_area(R,R,r) + 2*sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b8":
        return (cov(2*R, 2*R, 2*R, R, R, R)+Qpi*r**2) / tile_area
    elif bin_name == "b9":
        return (QQ(.5)*R**2+3*Qpi*r**2) / (sqrt(3)*R**2)
    else:
        print("I don't know how to compute it...")
        return 0
    
def almost_compact_density_r_greater(bin_name, r, R):
    # !!! need upper bound on r however
    if bin_name == "b1":
        a = (2*Qpi-2*tight_angle((R,r,R))-2*tight_angle((r,r,R)))/2
        if RIF(r/(r+R)) > RIF(cos(a)): # r is too big
            
            return 0               
        return Qpi*(R**2 + r**2) / (2*d_area(R,R,r) + d_area(R,r,r) + sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b2":
        a = (2*Qpi - tight_angle((r,r,r)) - tight_angle((R,r,R)) - 2*tight_angle((r,r,R))) / 2
        if RIF(r/(r+R)) > RIF(cos(a)) : # r is too big
            return 0
        return Qpi*(R**2 + 2*r**2) / (4*d_area(R,R,r) + 2*sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b3":
        a = (Qpi - 2*tight_angle((r,R,r))) / 2
        if RIF(R/(r+R)) > RIF(cos(a)) : # r is too big
            return 0
        return Qpi*(2*r**2 + R**2) / (4*d_area(r,R,r) + 2*sin(a)*cos(a)*(r+R)**2) 
    elif bin_name == "b4":
        a = (2*Qpi-3*tight_angle((R,r,R)))/2
        if RIF(r/(r+R)) > RIF(cos(a)): # r is too big
            return 0
        return Qpi*(R**2 + r**2) / (3*d_area(R,R,r) + sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b5":
        a = (2*Qpi-2*tight_angle((r,r,r))-2*tight_angle((R,r,r)))/2
        if RIF(r/(r+R)) > RIF(cos(a)): # r is too big
            return 0
        return Qpi*(2*R**2 + 7*r**2) / (6*d_area(r,r,R) +6*d_area(r,r,r) + 6*sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b6":
        RR = r/bin_dict['b6'][0]
        a = arccos(r/(RR+r))
        return Qpi*(3*r**2 + QQ(.5)*R**2) / (sqrt(3)*(sin(a)*(r+RR))**2)
    elif bin_name == "b7":
        a = (2*Qpi-tight_angle((R,r,R))-2*tight_angle((r,r,R)))/2
        if RIF(r/(r+R)) > RIF(cos(a)): # r is too big
            return 0
        return Qpi*(QQ(.5)*R**2 + r**2) / (d_area(R,R,r) + d_area(r,r,R) + sin(a)*cos(a)*(r+R)**2) 
    elif bin_name == "b8":
        a = Qpi-tight_angle((R,r,R))
        if RIF((2*r)/(r+R)) > RIF(cos(a)): # r is too big, or too close to b1
            return 0
        return Qpi*(QQ(.5)*R**2 + r**2) / (2*d_area(R,R,r) + sin(a)*cos(a)*(r+R)**2)
    elif bin_name == "b9":
        a = (2*Qpi-2*tight_angle((R,r,r))-tight_angle((r,r,r)))/2
        if  RIF(r/(r+R)) > RIF(cos(a)) : # r is too big
            return 0
        return Qpi*(QQ(.5)*R**2 + 3*r**2) / ( sqrt(3) * (sin(a)*(r+R))**2 )
    else:
        print("I don't know how to compute it...")
        return 0
    
        
def classify(bin_name,case, r, R, d):
    bin_r = bin_dict[bin_name][0]
    bin_d = bin_dict[bin_name][1]
    sign = ''
    if  r/R <= bin_r:
        sign = "r <= binr"
        d_new = almost_compact_density_r_smaller(bin_name, r, R)
    else:
        sign = "r > binr"
        d_new = almost_compact_density_r_greater(bin_name, r, R)

    if d_new - bin_d > QQ(0.00001):
        print(str(case)+"  Non compact bin > compact ???? something is wrong"+str((r,R,d,bin_name, (RIF(d_new).endpoints(), RIF(bin_dict[bin_name][1]).endpoints()))))
        return('','',0)
    if d_new > d:
        return("counter example",sign, d_new)
    else:
        return("no",sign, d_new)

 


                       
#########################################################################################
# Main
#########################################################################################

ternary_init() # it's long...
print("Ternary init done")

candidates = {}
for name in bin_names:
    candidates.update({name: ternary_with_similar_ratio(name)})

results_dict = {}
found_cases = set()
    
# consists of bin_name:[(ter_case, (ter_d, d_new), rtype, r, R, sign)]
for bin_name in bin_names:
    results = []
    print("treating "+bin_name+" ...")
    for rtype,candidates_list in candidates[bin_name]:
        for candidate in candidates_list:
            case = candidate['case']
            ter_d = candidate['ter_d']
            r,R=candidate['r'],candidate['R']
            candidate_type, sign,d_new = classify(bin_name,case,r,R,ter_d)
            if candidate_type == "counter example":
                results.append((case, (ter_d, d_new), rtype, r, R, sign))
                found_cases.add(case)
    results.sort(key = lambda x: x[0])
    results_dict.update({bin_name:results})


found_list = list(found_cases)
found_list.sort()
print("Found cases: "+str(found_list))

found_cases_info = {}
for case in found_list:
    l = []
    for bin_name in bin_names:
        for (ter_case, fdensities, rtype,r,R,sign) in results_dict[bin_name]:
            if ter_case==case:
                l+=[{'bin_name':bin_name, 'ter_d':fdensities[0], 'd_new':fdensities[1], 'rtype':rtype, 'r':r,'R':R,'sign':sign[2]}]
                break
    found_cases_info.update({case:l})


best_bin_for_found = {}

for case in found_cases_info.keys():
    l = found_cases_info[case]
    l.sort(key = lambda x: n(x['d_new']))
    binn=l[-1]
    best_bin_for_found.update({case:binn})


ternary_for_bin = {binname:[] for binname in bin_names}

not_saturated = [24,28,29,30,31,32,33,37,38,39,40,41,42,43,44]
found_saturated = [x for x in found_list if not x in not_saturated]

for case in found_saturated:
    bin_name=best_bin_for_found[case]['bin_name']
    l = ternary_for_bin[bin_name]
    ternary_for_bin.update({bin_name:(l+[case])})


bin_names = [x for x in bin_names if len(ternary_for_bin[x])>0]

for b in ternary_for_bin.keys():
    print(b+": "+str(ternary_for_bin[b]))
