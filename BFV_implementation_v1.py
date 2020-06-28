import numpy as np
import math
from scipy import stats
import random
random.seed(9823)

############################################

#discrete Gaussian supported on [-10*sigma, 10*sigma]
def disc_gaus(sigma):
    B = round(10*sigma)
    x = np.arange(-B,B+1)
    pmf = np.exp(-math.pi*x*x/sigma/sigma)
    pmf = pmf/np.sum(pmf)
    return stats.rv_discrete(values=(x, pmf))

#polynomial ring multiplication modulo x^d+1
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

#LPR.ES scheme based on RLWE problem
def key_gen():
    s = np.random.randint(2, size=d)
    a = np.random.randint(low=-q/2+1, high=q/2+1, size=d)
    e = chi.rvs(size=d)
    p0 = cent_red(-ring_mul(a,s)-e, q)
    p1 = a
    return s, [p0,p1]

def encrypt(pk, m):
    u = np.random.randint(2, size=d)
    e1 = chi.rvs(size=d)
    e2 = chi.rvs(size=d)
    c0 = cent_red(ring_mul(pk[0],u)+e1+delta*m, q)
    c1 = cent_red(ring_mul(pk[1],u)+e2, q)
    return [c0,c1]

def decrypt(sk, ct):
    c0 = ct[0]
    c1 = ct[1]
    return cent_red(nearest_round(cent_red(c0+ring_mul(c1,sk), q)*t/q), t)

#homomorphic addition
def add(ct1, ct2):
    ct_sum_0 = cent_red(ct1[0]+ct2[0], q)
    ct_sum_1 = cent_red(ct1[1]+ct2[1], q)
    return [ct_sum_0, ct_sum_1]

#homomorphic multiplication
#with ct1 and ct2 two ciphertexts of two messages, it returns the ciphertexts for the ring product of the two messages
def mul(ct1, ct2, v):
    c0 = cent_red(nearest_round(ring_mul(ct1[0],ct2[0])*t/q), q)
    c1 = cent_red(nearest_round((ring_mul(ct1[0],ct2[1])+ring_mul(ct1[1],ct2[0]))*t/q), q)
    c2 = cent_red(nearest_round(ring_mul(ct1[1],ct2[1])*t/q), q)
    if v == 1:
        return relin_1(c0,c1,c2)
    elif v == 2:
        return relin_2(c0,c1,c2)
    else:
        print('Wrong relinearisation version.')

#relinearisation version 1
def relin_1_keygen():
    rlk_v1_0 = []
    rlk_v1_1 = []
    for i in range(l+1):
        ai = np.random.randint(low=-q/2+1, high=q/2+1, size=d)
        ei = chi.rvs(size=d)
        rlk_v1_0.append(cent_red(-ring_mul(ai,sk)-ei+T**i*ring_mul(sk,sk), q))
        rlk_v1_1.append(ai)
    return [rlk_v1_0, rlk_v1_1]

def relin_1(c0,c1,c2):
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
        sum_0 = sum_0+ring_mul(rlk_v1[0][i],c2i)
        sum_1 = sum_1+ring_mul(rlk_v1[1][i],c2i)
    return [cent_red(c0+sum_0, q), cent_red(c1+sum_1, q)]

#relinearisation version 2
def relin_2_keygen():
    a = np.random.randint(low=-p*q/2+1, high=p*q/2+1, size=d)
    e = chi_prime.rvs(size=d)
    return [cent_red(-ring_mul(a,sk)-e+p*ring_mul(sk,sk), p*q), a]

def relin_2(c0,c1,c2):
    c2_0 = cent_red(nearest_round(ring_mul(c2,rlk_v2[0])/p), q)
    c2_1 = cent_red(nearest_round(ring_mul(c2,rlk_v2[1])/p), q)
    return [cent_red(c0+c2_0,q), cent_red(c1+c2_1,q)]

##################################################################

#initialisation
d = 8
q = 3**12
t = 10
delta = math.floor(q/t)
sk, pk = key_gen()

sigma = 1
chi = disc_gaus(sigma)

T = 3
l = math.floor(math.log(q,T))
rlk_v1 = relin_1_keygen()

sigma_prime = 10
chi_prime = disc_gaus(sigma_prime)
p = 1000
rlk_v2 = relin_2_keygen()

###########################################################################

#test1 
m = np.array([1,2,3,4,5,4,3,2])
ct = encrypt(pk, m)
print('ciphertext:')
print(ct)
print('decrypted:')
print(decrypt(sk, ct))

#test2 addition
m1 = np.array([1,2,3,4,5,4,3,2])
m2 = np.array([1,2,3,4,5,4,3,2])
ct1 = encrypt(pk,m1)
ct2 = encrypt(pk,m2)
print('true sum:')
print(cent_red(m1+m2,t))
print('decrypted:')
print(decrypt(sk, add(ct1,ct2)))

#test3 multiplication
m1 = np.array([1,0,1,4,5,-4,5,2])
m2 = np.array([1,2,-1,2,5,4,-4,2])
ct1 = encrypt(pk,m1)
ct2 = encrypt(pk,m2)
print('true m1*m2:')
print(cent_red(ring_mul(m1,m2), t))
print('decrypted:')
print(decrypt(sk, mul(ct1,ct2,2)))

#test4 more layers of multiplication
m1 = np.array([1,-2,1,1,5,-4,2,2])
m2 = np.array([0,2,-1,-4,5,4,-4,0])
m3 = np.array([1,1,1,4,2,-2,5,2])
m4 = np.array([1,1,-4,3,5,1,-1,1])
ct1 = encrypt(pk,m1)
ct2 = encrypt(pk,m2)
ct3 = encrypt(pk,m3)
ct4 = encrypt(pk,m4)
m1xm2 = cent_red(ring_mul(m1,m2), t)
m3xm4 = cent_red(ring_mul(m3,m4), t)
print('true m1*m2*m3*m4')
print(cent_red(ring_mul(m1xm2,m3xm4), t))
print('decrypted')
print(decrypt(sk, mul(mul(ct1,ct2,1),mul(ct3,ct4,1),1)))
