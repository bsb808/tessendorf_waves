
from math import *

def dispersion(n_prime, m_prime):
    N = 10
    
    w_0 = 2.0 * pi / 200.0;
    kx = pi * (2 * n_prime - N) / length;
    kz = M_PI * (2 * m_prime - N) / length;
    return floor(sqrt(g * sqrt(kx * kx + kz * kz)) / w_0) * w_0;

def phillips(n_prime, m_prime):
    vector2 k(M_PI * (2 * n_prime - N) / length,
          M_PI * (2 * m_prime - N) / length);
    float k_length  = k.length();
    if (k_length < 0.000001) return 0.0;
 
    float k_length2 = k_length  * k_length;
    float k_length4 = k_length2 * k_length2;
 
    float k_dot_w   = k.unit() * w.unit();
    float k_dot_w2  = k_dot_w * k_dot_w;
 
    float w_length  = w.length();
    float L         = w_length * w_length / g;
    float L2        = L * L;
     
    float damping   = 0.001;
    float l2        = L2 * damping * damping;
 
    return A * exp(-1.0f / (k_length2 * L2)) / k_length4 * k_dot_w2 * exp(-k_length2 * l2);
