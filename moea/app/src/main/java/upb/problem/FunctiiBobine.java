package upb.problem;

import static java.lang.Math.*;
import static upb.problem.MatlabFunctions.ellipke;
import static upb.problem.MatlabFunctions.linspace;

class FunctiiBobine {

    /**
     * Calculeaza intensitatea campului magnetic produs de o spira filiforma circulara
     * de raza a [m], parcursa de curentul i [A] aflata in planul de cota z0 (unde axa
     * Oz a unui sistem de coordonate cilindric este perpendiculara pe planul spirei si
     * trece prin centrul spirei) Campul se calculeaza in punctul de coordonate coord
     * (intr-un sistem de coordonate cilindric; unghiul nu conteaza, problema are o
     * simetrie cilindrica. Campul are doar 2 componente, dupa r si z
     * <p>
     * Data ultimei modificari: 23 februarie 2023
     *
     * @param a     = raza spirei [m]
     * @param i     = intensitatea curentului [A]
     * @param z0    = cota planului in care se afla spira [m]
     * @param coord = este o matrice cu mai multe linii (cel putin 1) si 2 coloane;
     *              fiecare linie corespunde unui punct
     *              coord(k,1) reprezinta raza r [m]
     *              coord(k,2) reprezinta cota z [m]
     * @param H     = componentele intensitatii campului magnetic [A/m]
     *              este o matrice cu un numar de linii egal cu cel al vectorului coord
     *              si un numar de coloane egal cu 2
     *              H(k,1) reprezinta Hr
     *              H(k,2) reprezinta Hz
     * @author Gabriela Ciuprina
     */
    static void campHoriunde_spira(double a, double i, double z0, double[][] coord, double tol, double[][] H) {
        double a2 = pow(a, 2);
        double C = i / PI;
        for (int j = 0; j < coord.length; j++) {
            double r = coord[j][0];
            double z = coord[j][1];
            double r2 = pow(r, 2);
            double zmz0 = z - z0;
            double z2 = pow(zmz0, 2);
            double doiar = 2 * a * r;

            double a2pr2pz2 = a2 + r2 + z2;
            double alpha2 = a2pr2pz2 - doiar;
            double beta2 = a2pr2pz2 + doiar;
            double doialpha2beta = 2 * alpha2 * sqrt(beta2);
            double k2 = 1 - alpha2 / beta2;

            double[] KE = ellipke(k2, tol);
            double a2mr2mz2 = a2 - r2 - z2;
            double alpha2K = alpha2 * KE[0];
            double Hr = C * zmz0 / doialpha2beta / r * (a2pr2pz2 * KE[1] - alpha2K);
            double Hz = C / doialpha2beta * (a2mr2mz2 * KE[1] + alpha2K);
            H[j][0] += Hr;
            H[j][1] += Hz;
        }
    }

    static void campHoriunde_spira(double a, double i, double z0, double z0_2, double[][] coord, double tol, double[][] H) {
        double a2 = pow(a, 2);
        double C = i / PI;
        for (int j = 0; j < coord.length; j++) {
            double r = coord[j][0];
            double z = coord[j][1];
            double r2 = pow(r, 2);
            double zmz0 = z - z0;
            double zmz0_2 = z - z0_2;
            double z2 = pow(zmz0, 2);
            double z2_2 = pow(zmz0_2, 2);
            double doiar = 2 * a * r;

            double a2pr2pz2 = a2 + r2 + z2;
            double a2pr2pz2_2 = a2 + r2 + z2_2;
            double alpha2 = a2pr2pz2 - doiar;
            double alpha2_2 = a2pr2pz2_2 - doiar;
            double beta2 = a2pr2pz2 + doiar;
            double beta2_2 = a2pr2pz2_2 + doiar;
            double doialpha2beta = 2 * alpha2 * sqrt(beta2);
            double doialpha2beta_2 = 2 * alpha2_2 * sqrt(beta2_2);
            double k2 = 1 - alpha2 / beta2;
            double k2_2 = 1 - alpha2_2 / beta2_2;

            double[] KE = ellipke(k2, tol);
            double[] KE_2 = ellipke(k2_2, tol);
            double a2mr2mz2 = a2 - r2 - z2;
            double a2mr2mz2_2 = a2 - r2 - z2_2;
            double alpha2K = alpha2 * KE[0];
            double alpha2K_2 = alpha2_2 * KE_2[0];
            double Hr = C * zmz0 / doialpha2beta / r * (a2pr2pz2 * KE[1] - alpha2K);
            double Hr_2 = C * zmz0_2 / doialpha2beta_2 / r * (a2pr2pz2_2 * KE_2[1] - alpha2K_2);
            double Hz = C / doialpha2beta * (a2mr2mz2 * KE[1] + alpha2K);
            double Hz_2 = C / doialpha2beta_2 * (a2mr2mz2_2 * KE_2[1] + alpha2K_2);
            H[j][0] += Hr + Hr_2;
            H[j][1] += Hz + Hz_2;
        }
    }

    /**
     * Calculeaza intensitatea campului magnetic produs de un solenoid infinit subtire
     * de raza a, cota minima zmin si lungime h, parcurs de o solenatie NI Calculul se
     * face prin superpozitie, se apeleaza functia care evalueaza integrale eliptice
     * Campul se calculeaza in punctul de coordonate coord (intr-un sistem de
     * coordonate cilindric; unghiul nu conteaza, problema are o simetrie cilindrica
     * Campul are doar 2 componente, dupa r si z
     * <p>
     * Data ultimei modificari: 23 februarie 2023
     *
     * @param a       raza spirei [m]
     * @param coord   este o matrice cu mai multe linii (cel putin 1) si 2 coloane;
     *                fiecare linie corespunde unui punct
     *                coord(k,1) reprezinta raza r [m]
     *                coord(k,2) reprezinta cota z [m]
     * @param nrspire numarul de spire in care se discretizeaza solenoidul
     * @param tol     toleranta pentru integralele eliptice
     * @param H       = componentele intensitatii campului magnetic [A/m]
     *                este o matrice cu un numar de linii egal cu cel al vectorului coord
     *                si un numar de coloane egal cu 2
     *                H(k,1) reprezinta Hr
     *                H(k,2) reprezinta Hz
     * @author Gabriela Ciuprina
     * Data ultimei modificari: 23 februarie 2023
     */
    static void campHoriunde_solenoidInfinitSubtire_superpozitie(
            double a,
            double zmin,
            double h,
            double NI,
            double[][] coord,
            int nrspire,
            double tol,
            double[][] H
    ) {
        double zmax = zmin + h;
        double[] vector_z0 = linspace(zmin, zmax, nrspire);
        double crtspira = NI / nrspire;
        for (int k = 0; k < nrspire; k++) {
            campHoriunde_spira(a, crtspira, vector_z0[k], coord, tol, H);
        }
    }

    static void campHoriunde_solenoidInfinitSubtire_superpozitie(
            double a,
            double zmin,
            double zmin2,
            double h,
            double NI,
            double[][] coord,
            int nrspire,
            double tol,
            double[][] H
    ) {
        double zmax = zmin + h;
        double zmax2 = zmin2 + h;
        double[] vector_z0 = linspace(zmin, zmax, nrspire);
        double[] vector_z0_2 = linspace(zmin2, zmax2, nrspire);
        double crtspira = NI / nrspire;
        for (int k = 0; k < nrspire; k++) {
            campHoriunde_spira(a, crtspira, vector_z0[k], vector_z0_2[k], coord, tol, H);
        }
    }

    /**
     * Calculeaza intensitatea campului magnetic produs de o bobina de raza interioara
     * a1[ m], si raza exterioara a2 [m], parcursa de solenatia NI [A*spire], avand
     * cota minima zmin (unde axa Oz a unui sistem de coordonate cilindric este coaxial
     * cu bobina) Campul se calculeaza in punctul de coordonate coord (intr-un sistem
     * de coordonate cilindric; unghiul nu conteaza, problema are o simetrie cilindrica
     * Campul are doar 2 componente, dupa r si z
     * <p>
     * Data ultimei modificari: 23 februarie 2023
     *
     * @param a1                 raza interioara a bobinei [m] (raza cilindrului interior)
     * @param a2                 raza exterioara a bobinei [m]
     * @param zmin               cota minima a solenoidului [m]
     * @param h                  lungimea solenoidului [m]
     * @param NI                 solenatia = intensitatea curentului [A] * numarul de spire
     * @param coord              este o matrice cu mai multe linii (cel putin 1) si 2 coloane;
     *                           fiecare linie corespunde unui punct
     *                           coord(k,1) reprezinta raza r [m]
     *                           coord(k,2) reprezinta cota z [m]
     * @param noSInfSub          numarul de solenoizi infinit subtiri pentru discretizare
     * @param nrspirePerSolenoid numarul de spire in care se discretizeaza fiecare solenoidul
     * @param tol                toleranta pentru integralele eliptice
     * @return H                  = componentele intensitatii campului magnetic [A/m]
     * este o matrice cu un numar de linii egal cu cel al vectorului coord
     * si un numar de coloane egal cu 2
     * H(k,1) reprezinta Hr
     * H(k,2) reprezinta Hz
     * @author Gabriela Ciuprina
     */
    static double[][] campHoriunde_bobina_supSpElem(
            double a1,
            double a2,
            double zmin,
            double h,
            double NI,
            double[][] coord,
            int noSInfSub,
            int nrspirePerSolenoid,
            double tol
    ) {
        double[][] H = new double[coord.length][2];
        double NIsol = NI / noSInfSub;
        double[] a_vector = linspace(a1, a2, noSInfSub);
        for (int k = 0; k < noSInfSub; k++) {
            campHoriunde_solenoidInfinitSubtire_superpozitie(a_vector[k], zmin, h, NIsol, coord, nrspirePerSolenoid, tol, H);
        }
        return H;
    }

    static double[][] campHoriunde_bobina_supSpElem(
            double a1,
            double a2,
            double zmin,
            double zmin2,
            double h,
            double NI,
            double[][] coord,
            int noSInfSub,
            int nrspirePerSolenoid,
            double tol
    ) {
        double[][] H = new double[coord.length][2];
        double NIsol = NI / noSInfSub;
        double[] a_vector = linspace(a1, a2, noSInfSub);
        for (int k = 0; k < noSInfSub; k++) {
            campHoriunde_solenoidInfinitSubtire_superpozitie(a_vector[k], zmin, zmin2, h, NIsol, coord, nrspirePerSolenoid, tol, H);
        }
        return H;
    }

    /**
     * Calculeaza (cu formula analitica simpla) intensitatea campului magnetic produs
     * de un solenoid infinit subtire, de raza a [m], parcurs de solenatia NI
     * [A*spire], avand cota minima zmin (unde axa Oz a unui sistem de coordonate
     * cilindric este coaxial cu solenoidul) Campul se calculeaza in punctul de cota z
     * (poate fi vector). Pe axa campul are doar componenta orientata dupa z
     * <p>
     * Data ultimei modificari: 23 februarie 2023
     *
     * @param a    raza spirei [m]
     * @param NI   solenatia = intensitatea curentului [A] * numarul de spire
     * @param zmin cota minima a solenoidului [m]
     * @param h    lungimea solenoidului [m]
     * @param z    cota punctului de pe axa unde se calculeaza campul [m] (poate fi un vector)
     * @param Hz   = componenta dupa z a intensitatii campului magnetic [A/m] (daca z este un vector, atunci acest Hz este un vector)
     * @author Gabriela Ciuprina
     */
    static void campHaxa_solenoidInfinitSubtire(double a, double zmin, double h, double NI, double[] z, double[] Hz) {
        double a2 = pow(a, 2);
        for (int i = 0; i < z.length; i++) {
            double zmzmin = z[i] - zmin;
            double cos_alpha1 = zmzmin / sqrt(pow(zmzmin, 2) + a2);
            double zmax = zmin + h;
            double zmaxmz = zmax - z[i];
            double cos_alpha2 = zmaxmz / sqrt(pow(zmaxmz, 2) + a2);
            Hz[i] += NI / (2 * h) * (cos_alpha1 + cos_alpha2);
        }
    }

    static void campHaxa_solenoidInfinitSubtire(double a, double zmin, double zmin2, double h, double NI, double[] z, double[] Hz) {
        double a2 = pow(a, 2);
        double NIover2h = NI / (2 * h);
        for (int i = 0; i < z.length; i++) {
            double zmzmin = z[i] - zmin;
            double zmzmin2 = z[i] - zmin2;
            double cos_alpha1 = zmzmin / sqrt(pow(zmzmin, 2) + a2);
            double cos_alpha1_2 = zmzmin2 / sqrt(pow(zmzmin2, 2) + a2);
            double zmax = zmin + h;
            double zmax2 = zmin2 + h;
            double zmaxmz = zmax - z[i];
            double zmaxmz2 = zmax2 - z[i];
            double cos_alpha2 = zmaxmz / sqrt(pow(zmaxmz, 2) + a2);
            double cos_alpha2_2 = zmaxmz2 / sqrt(pow(zmaxmz2, 2) + a2);

            Hz[i] += NIover2h * (cos_alpha1 + cos_alpha2 + cos_alpha1_2 + cos_alpha2_2);
        }
    }

    /**
     * Calculeaza intensitatea campului magnetic produs de o bobina de raza interioara
     * a1[m], si raza exterioara a2 [m], parcursa de solenatia NI [A*spire], avand cota
     * minima zmin (unde axa Oz a unui sistem de coordonate cilindric este coaxial
     * cu bobina)
     * Campul se calculeaza in punctul de cota z (poate fi vector). Pe axa campul are
     * doar componenta orientata dupa z Calculul campului se va face prin superpozitie
     * de campuri produse de solenoizi infinit subtiri, campurile fiind calculate cu
     * formula analitica simpla
     * <p>
     * Data ultimei modificari: 23 februarie 2023
     *
     * @param a1        raza interioara a bobinei [m] (raza cilindrului interior)
     * @param a2        raza exterioara a bobinei [m]
     * @param NI        solenatia = intensitatea curentului [A] * numarul de spire
     * @param zmin      cota minima a solenoidului [m]
     * @param h         lungimea solenoidului [m]
     * @param z         cota punctului de pe axa unde se calculeaza campul [m]
     *                  (poate fi un vector)
     * @param noSInfSub numarul de solenoizi infinit subtiri in care se va
     *                  discretiza bobina;
     * @return Hz        = componenta dupa z a intensitatii campului magnetic [A/m]
     * (daca z este un vector, atunci acest Hz este un vector)
     * @author Gabriela Ciuprina
     */
    static double[] campHaxa_bobina_supSInfSub(double a1, double a2, double zmin, double h, double NI, double[] z, int noSInfSub) {
        double[] Hz = new double[z.length];
        double NIsol = NI / noSInfSub;
        double[] a_vector = linspace(a1, a2, noSInfSub);
        for (int k = 0; k < noSInfSub; k++) {
            campHaxa_solenoidInfinitSubtire(a_vector[k], zmin, h, NIsol, z, Hz);
        }
        return Hz;
    }

    static double[] campHaxa_bobina_supSInfSub(double a1, double a2, double zmin, double zmin2, double h, double NI, double[] z, int noSInfSub) {
        double[] Hz = new double[z.length];
        double NIsol = NI / noSInfSub;
        double[] a_vector = linspace(a1, a2, noSInfSub);
        for (int k = 0; k < noSInfSub; k++) {
            campHaxa_solenoidInfinitSubtire(a_vector[k], zmin, zmin2, h, NIsol, z, Hz);
        }
        return Hz;
    }

}