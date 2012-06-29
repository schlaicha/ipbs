#include <gsl/gsl_sf_ellint.h>

template<class V, class W> 
W E_field_cartesian(V r, V r_prime) {
//    std::cout << "E Field cartesian!" << std::endl;
              V dist = r - r_prime;
              W E_field;
              E_field = dist;
#if     GRIDDIM == 2
                E_field /= dist.two_norm() * dist.two_norm();
                E_field *= 2;
#elif   GRIDDIM == 3
                E_field /= dist.two_norm() * dist.two_norm() * dist.two_norm();
#else
                DUNE_THROW(Dune::NotImplemented, "IPBS only knows about dim=2 and dim=3");
#endif                  
                return E_field;
}

template<class V, class W>
W E_field_cylinder(V r, V r_prime) {
//    std::cout << "E Field cylinder!" << std::endl;

              V dist = r - r_prime;
                double a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
                double b = 4.0 * r[1] * r_prime[1];
                double k = sqrt(b/a);
                double E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
                double K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
                W E_field(0.);
                E_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
                E_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1])
                               / ((a-b)*sqrt(a)*b) * E
                             + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
                E_field*=-2;
                return E_field;
}

template<class V, class W>
W E_field_cylinder_with_mirror(V r, V r_prime) {
    V r_prime2 = r_prime;
    r_prime2[0] *= -1;

    //std::cout << "r " << r << " r' " << r_prime << " r2' " << r_prime2 << " E_field: " <<  E_field_cylinder<V, W>(r, r_prime)  << " mirrored: " <<  E_field_cylinder<V, W>(r, r_prime2) << std::endl;
    return ( E_field_cylinder<V, W>(r, r_prime)  + E_field_cylinder<V, W>(r, r_prime2) );
}

template< class V, class W> 
W E_field(V r, V r_prime, int geometry) 
{
    if (geometry == 0) {
        return E_field_cartesian<V, W> (r, r_prime);
    } else if (geometry == 1) {
        return E_field_cylinder<V, W> (r, r_prime);
    } else if (geometry == 2) {
        return E_field_cylinder_with_mirror<V, W>(r, r_prime);
    } else {
        return W(0);
    }
}
