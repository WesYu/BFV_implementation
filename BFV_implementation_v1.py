import numpy as np
import math
from scipy import stats
import random
random.seed(624874)

###############################################

#discrete Gaussian supported on [-10*sigma, 10*sigma]
def disc_gaus(sigma):
    B = round(10*sigma)
    x = np.arange(-B,B+1)
    pmf = np.exp(-math.pi*x*x/sigma/sigma)
    pmf = pmf/np.sum(pmf)
    return stats.rv_discrete(values=(x, pmf))

#polynomial ring multiplication
def ring_mul(a1, a2):
    d = a1.shape[0]
    result = np.zeros(d)
    for i in range(d):
        for j in range(d):
            if i+j>=d:
                result[i+j-d] = result[i+j-d]-a1[i]*a2[j]
            else:
                result[i+j] = result[i+j]+a1[i]*a2[j]
    return result

#move coefficients to (-q/2, q/2]
def cent_red(a, q):
    d = a.shape[0]
    result = np.zeros(d)
    for i in range(d):
        result[i] = a[i]+q*math.floor(1/2-a[i]/q)
    return result

#round the coefficients to the nearest integer
def nearest_round(a):
    d = a.shape[0]
    result = np.zeros(d)
    for i in range(d):
        result[i] = int(round(a[i]))
    return result

###############################################

#LPR.ES scheme based on RLWE problem
#optimised by sampling s and u from R2
d = 16
q = 10000
t = 3
delta = math.floor(q/t)

sigma = 1
chi = disc_gaus(sigma)

#secrete key generation
s = np.random.randint(2, size=d)

#public key generation
a = np.random.randint(low=-q/2+1, high=q/2+1, size=d)
e = chi.rvs(size=d)
p0 = cent_red(-ring_mul(a,s)-e, q)
p1 = a

#encryption
m = np.random.randint(low=-t/2+1, high=t/2+1, size=d)
u = np.random.randint(2, size=d)
e1 = chi.rvs(size=d)
e2 = chi.rvs(size=d)
c0 = cent_red(ring_mul(p0,u)+e1+delta*m, q)
c1 = cent_red(ring_mul(p1,u)+e2, q)

#decryption
decrypt = cent_red(nearest_round(cent_red(c0+ring_mul(c1,s), q)*t/q), t)
print(m-decrypt)

###############################################

#homomorphic addition
#encryption
m1 = np.random.randint(low=-t/2+1, high=t/2+1, size=d)
m2 = np.random.randint(low=-t/2+1, high=t/2+1, size=d)
c10 = cent_red(ring_mul(p0,u)+e1+delta*m1, q)
c11 = cent_red(ring_mul(p1,u)+e2, q)
c20 = cent_red(ring_mul(p0,u)+e1+delta*m2, q)
c21 = cent_red(ring_mul(p1,u)+e2, q)

#homomorphic addition
add_c0 = cent_red(c10+c20, q)
add_c1 = cent_red(c11+c21, q)

#decryption
add_decrypt = cent_red(nearest_round(cent_red(add_c0+ring_mul(add_c1,s), q)*t/q), t)

#test correctness
plain_add = cent_red(m1+m2, t)
print(plain_add-add_decrypt)

############################################

#homomorphic multiplication with relinearisation version 1
#relinearisation key generation
T = 10
l = math.floor(math.log(q,T))
rlk0_v1 = []
rlk1_v1 = []
for i in range(l+1):
    ai = np.random.randint(low=-p*q/2+1, high=p*q/2+1, size=d)
    ei = chi.rvs(size=d)
    rlk0_v1.append(cent_red(-ring_mul(ai,s)-ei+T**i*ring_mul(s,s), q))
    rlk1_v1.append(ai)
    
#homomorphic multiplication
c0 = cent_red(nearest_round(ring_mul(c10,c20)*t/q), q)
c1 = cent_red(nearest_round((ring_mul(c10,c21)+ring_mul(c11,c20))*t/q), q)
c2 = cent_red(nearest_round(ring_mul(c11,c21)*t/q), q)
c2_temp = c2
c2_T = []
sum_0 = np.array(d)
sum_1 = np.array(d)
for i in range(l+1):
    c2i = np.zeros(d)
    for j in range(d):
        if c2_temp[j]>=0:
            c2i[j] = c2_temp[j]%T
        else:
            c2i[j] = -(-c2_temp[j])%T
        c2_temp[j] = (c2_temp[j]-c2i[j])/T
    c2_T.append(c2i)
    sum_0 = sum_0+ring_mul(rlk0_v1[i],c2i)
    sum_1 = sum_1+ring_mul(rlk1_v1[i],c2i)
mul_c0_v1 = cent_red(c0+sum_0, q)
mul_c1_v1 = cent_red(c1+sum_1, q)

#decryption
mul_decrypt_v1 = cent_red(nearest_round(cent_red(mul_c0_v1+ring_mul(mul_c1_v1,s), q)*t/q), t)

#test accuracy
plain_mul = cent_red(ring_mul(m1,m2), t)
print(plain_mul-mul_decrypt_v1)

#########################################################

#homomorphic multiplication with relinearisation version 2
#chi-prime distribution
sigma = 100
chi_prime = disc_gaus(sigma)

#relinearisation key generation
p = 30
a = np.random.randint(low=-p*q/2+1, high=p*q/2+1, size=d)
e = chi_prime.rvs(size=d)
rlk0_v2 = cent_red(-ring_mul(a,s)-e+p*ring_mul(s,s), p*q)
rlk1_v2 = a

#homomorphic multiplication
c2_0 = cent_red(nearest_round(ring_mul(c2,rlk0_v2)/p), q)
c2_1 = cent_red(nearest_round(ring_mul(c2,rlk1_v2)/p), q)
mul_c0_v2 = cent_red(c0+c2_0, q)
mul_c1_v2 = cent_red(c1+c2_1, q)

#decryption
mul_decrypt_v2 = cent_red(nearest_round(cent_red(mul_c0_v2+ring_mul(mul_c1_v2,s), q)*t/q), t)

#test accuracy
print(plain_mul-mul_decrypt_v2)
