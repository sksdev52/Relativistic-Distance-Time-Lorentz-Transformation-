import math
import scipy
import scipy.constants as const
Pi = math.pi

light_speed = const.physical_constants['speed of light in vacuum'] # Returns a tuple whose first entry is the speed of light in metres per second
#print("Speed in light in vacuum", light_speed[0])
ltsp = light_speed[0]      # Not using c for speed of light to avoid any chance of confusing it with any other loop variable c
tvl = 0.999 * ltsp # travel velocity

def get_beta_gamma (velocity, receding = True):
    v = velocity
    to_send_beta_gamma = []
    if v ==0: # If velocity is 0
        return [1,1] # beta and gamma remain 1
    elif receding == False:
        v = -v
    beta = v/ltsp
    to_send_beta_gamma.append(beta)
    gamma = 1/math.sqrt(1-(beta**2))
    to_send_beta_gamma.append(gamma)
    return to_send_beta_gamma

def relativistic_beta_gamma_addition( vlist = [], rec = []): # List of relative velocities and list to indicate if the primed observer is receding or approaching
    '''If two bodies are approaching each other at relative velocities 'receding' should be set to True'''
    to_send_beta_gamma = []
    if vlist[0] == 0:
        beta_initial = 0
    elif rec[0] == False:
        vlist[0] = -vlist[0]
    else:
        beta_initial = vlist[0]/ltsp
    #print("Initial beta", beta_initial) # For testing purpose
    for i in range (1, len(vlist)):
            if vlist[i]==0: # If some intermediatory velocity is 0
                beta_final = beta_initial # Set the final beta to initial beta
            elif rec[i] == False:
                vlist[i] = -vlist[i]
                
            bnum = beta_initial + (vlist[i]/ltsp)# beta numerator
            if bnum == 0:
                beta_final = beta_initial
                #print(f"Beta Loop {i}", beta_final) # For testing
            else:
                bdm = 1 + (beta_initial * (vlist[i]/ltsp))# beta denominator
                beta_final = bnum/bdm
                #print(f"Beta Loop {i}", beta_final) # For testing
            beta_initial = beta_final
    if beta_final == 0:
        beta_final = 1 # Prevent division by 0 in a future function call.
    gamma_final = 1/math.sqrt(1-(beta_final**2))
    if gamma_final <= 1.0000001: # In case gamma final is not exactly 1 but with a very low decimal value
        gamma_final = 1.0
    to_send_beta_gamma.append(beta_final)
    to_send_beta_gamma.append(gamma_final)
    
    return to_send_beta_gamma
   


def vector_dot_vector(a=[], b=[]):
    if len(a) != len(b):
        print ("Vector Elements not equal to Matrix Columns")
        return
    else:
        sum = 0.0
        for e in range (len(a)): #iterate through each element
                sum += a[e] * b[e]   
        return sum

def lnz_back_transform_matrix_2D(beta_gamma = []): #arguments = gamma and beta, Lambda super mu sub nu, spcr = spherical co-ordinate angles in radians in a list, 1st value = polar angle, 2nd value =azimuthal angle
    # This transformation requires that one of the axis be aligned in the direction of motion
    bt = beta_gamma[0]    
    gm= beta_gamma[1]
    
    btm_to_send = [  # Forward transform matrix to send
                [gm, -(gm*bt)],
                [-(gm*bt), gm],
                 ]
    return btm_to_send



bg = get_beta_gamma (tvl, receding = True)
print("Beta, Gamma (Receding = True), along direction of motion: ", bg)
lft = lnz_back_transform_matrix_2D(bg) # Lorentz forward transform
uv_dir_m = [[1,0],[0,1]] #Unit vectors along direction of motion

trsuv= [] # transformed unit vectors distance
for e in range (len(lft)):
    element = vector_dot_vector(lft[e], uv_dir_m[e])
    trsuv.append(element)
    
print("Transformed unit vectors (time, distance) along direction of motion: ", trsuv)
print()
tm_change = trsuv[0]-1 # Change in the time vector, Final - initial . 1 can be changed to user desired value
ds_change = trsuv[1]-1 # Change in the distance vector, Final - initial . 1 can be changed to user desired value
spcr = [Pi/2, Pi/4] # Spherical co-ordinates, polar followed by azimuthal angles in radians. polar = pi/2 chosen to keep changes to 0 in direction i.e z direction

tm_contr = [] # Contribution to change by time components along different axis
contr = [] # Contribution to change by distance components along different axis

tm_ch_x = (tm_change * math.sin(spcr[0]) * math.sin(spcr[1])) # component contributed by x
tm_contr.append(tm_ch_x)
cmp_x= (ds_change * math.sin(spcr[0]) * math.sin(spcr[1])) # component contributed by x
contr.append(cmp_x)
print("Time, Distance change along x axis: ", tm_ch_x, cmp_x)
print()
tm_ch_y=(tm_change * math.sin(spcr[0]) * math.sin(spcr[1]))
tm_contr.append(tm_ch_y)
cmp_y= (ds_change * math.sin(spcr[0]) * math.sin(spcr[1])) # component contributed by y
contr.append(cmp_y)
print("Time, Distance  along y axis: ", tm_ch_y, cmp_y)
print()
tm_ch_z = (tm_change * math.cos(spcr[0]) )
tm_contr.append(tm_ch_z)
cmp_z= (ds_change * math.cos(spcr[0]) ) # component contributed by z
contr.append(cmp_z)
print("Time, Distance  along z axis: ", tm_ch_z, cmp_z)
print()

check_sum_tm = math.sqrt (tm_ch_x**2 + tm_ch_y**2 + tm_ch_z**2 )
print("Check to make sure that the time components add up to total change: ", check_sum_tm)
check_sum = math.sqrt (cmp_x**2 + cmp_y**2 + cmp_z**2 )
print("Check to make sure that the distance components add up to total change: ", check_sum)

final_unit_vectors_tm = []
final_tot_tm = 1 + tm_change
final_unit_vectors_tm.append(final_tot_tm) # # Total time might be needed. See final matrix below
final_unit_vectors_dist = []




for d in range (0, len(contr)):
    tm_cn = 1+ tm_contr[d] 
    final_unit_vectors_tm.append(tm_cn) # 1 can be replaced by the initial value. If change is negative it will be subtracted from the initial value
    lg_cn = 1+ contr[d] # Total distance not needed. Hence this list is 1 item less than time list
    final_unit_vectors_dist.append(lg_cn)

print()
print("Final time unit vectors", final_unit_vectors_tm)
print("Final distance unit vectors", final_unit_vectors_dist)

def lnz_back_transform_matrix_4D(new_tm_basis = [], new_dist_basis = []): # Vector components
    tm_tot = new_tm_basis[0]
    tm_x = new_tm_basis[1]
    tm_y = new_tm_basis[2]
    tm_z = new_tm_basis[3]
    
    ds_x = new_dist_basis[0]
    ds_y = new_dist_basis[1]
    ds_z = new_dist_basis[2]
    
    # Since the energy of a boosted particle is moved to the time vector, total time is required.
    matrix_4D_to_send = [           # This matrix needs to re-examined and modified if needed
        [tm_tot, ds_x, ds_y, ds_z], # Total time is used in this row
        [ds_x, tm_x, ds_y, ds_z],   # Components of time along each axis is used in subsequent rows
        [ds_y, ds_x, tm_y, ds_z],
        [ds_z, ds_x, ds_y, tm_z],
    ]
    return matrix_4D_to_send
print()
matrix4d = lnz_back_transform_matrix_4D(final_unit_vectors_tm , final_unit_vectors_dist)
print("Final 4D matrix")
for p in range (0, (len(matrix4d))):
    print (matrix4d[p])