import math
import scipy
import scipy.constants as const

light_speed = const.physical_constants['speed of light in vacuum'] # Returns a tuple whose first entry is the speed of light in metres per second
#print("Speed in light in vacuum", light_speed[0])
ltsp = light_speed[0]      # Not using c for speed of light to avoid any chance of confusing it with any other loop variable c
tvl = 0.5 * ltsp # travel velocity

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
        beta_final = 1 # Prevent multiplication by 0 in a future function call.
    gamma_final = 1/math.sqrt(1-(beta_final**2))
    if gamma_final <= 1.0000000: # In case beta is not exactly 0 but very very low value
        beta_final = 1.0
    to_send_beta_gamma.append(beta_final)
    to_send_beta_gamma.append(gamma_final)
    
    return to_send_beta_gamma
   

def event_primed_to_stationary_time_distance (vlist, rec, tm, dst): # vlist =relative velocity either one or a list of relavtive velocities
                                                                            #primed (moving) observer records an event at time tm  (ct) and at a distance (dst) 
    if vlist is list:
        btgm = relativistic_beta_gamma_addition( vlist , rec)
    else: 
       btgm = get_beta_gamma (vlist, rec)
    
    bt = btgm[0]
    gm = btgm[1]
    
    rtm = (tm + (gm * bt * dst))/gm
    rdst = (dst + (gm *bt * tm))/gm
    return rtm, rdst

def event_stationary_to_primed_time_distance (vlist, rec, tm, dst): # arguments: Time on the stationary clock, distance measured by stationary observer
    if vlist is list:
        btgm = relativistic_beta_gamma_addition( vlist , rec)
    else: 
        btgm = get_beta_gamma (vlist, rec)
    
    bt = btgm[0]
    gm = btgm[1]
    rtm = gm * (tm - (bt* dst))
    rdst = gm * (dst - (bt * tm))
    return rtm, rdst


bg = get_beta_gamma (tvl, receding = True)
print("Beta, Gamma (Receding = True)", bg)

bg2 = get_beta_gamma (tvl, receding = False)
print("Beta, Gamma (Receding = False)", bg2)
print


vel = [tvl, tvl, tvl, tvl, tvl, tvl, tvl] # Classic case of a spaceship which fires a missile at relativistic velocity which in turn fires a missile at relativistic velocity and so on
rec = [True,True,True,True,True, True, True]
ans = relativistic_beta_gamma_addition( vel, rec)
print ("All receding = True", ans)        
print()
vel2 = [tvl, tvl, tvl, tvl, tvl, tvl, tvl]
rec2 = [True,False,True,False, True, False, False]
ans2 = relativistic_beta_gamma_addition( vel2, rec2)
print ("Receding combination of true and False", ans2)