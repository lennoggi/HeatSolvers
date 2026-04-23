#include <cassert>
#include <cmath>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>

#include <hdf5.h>

#include "include/Declare_functions.hh"
#include "include/Check_parameters.hh"

#include "Parameters.hh"

using namespace std;


int main() {
    // Print some info
    #if (INITIAL_DATA == GAUSSIAN)
    cout << "Gaussian initial data" << endl
         << "Initial center: (" << X0     << ", " << Y0     << ")" << endl
         << "Initial sigma:  (" << SIGMAX << ", " << SIGMAY << ")" << endl
         << endl;
    #else
    #error "Invalid initial data"
    #endif


    // Build the initial data
    constexpr double dx = LX/static_cast<double>(NX);
    constexpr double dy = LY/static_cast<double>(NY);

    vector<double> x(NX), y(NY);  // NOTE: only used for output

    constexpr auto ntot = NX*NY;
    vector<double> u(ntot);

    for (auto i = decltype(NX){0}; i < NX; ++i) {
        x[i] = i*dx;
        const auto idx_i = i*NY;

        for (auto j = decltype(NY){0}; j < NY; ++j) {
            y[j] = j*dy;

            #if (INITIAL_DATA == GAUSSIAN)
            u[idx_i + j] = gaussian(x[i], y[j]);
            #else
            #error "Invalid initial data"
            #endif
        }
    }


    /* ------------------------
     * Initialize the HDF5 file
     * ------------------------ */
    const auto file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);


    // ***** Output frequency ****
    constexpr hsize_t one   = 1;
    const auto space_one_id = H5Screate_simple(1, &one, nullptr);
    assert(space_one_id >= 0);

    const auto out_freq_dset_id = H5Dcreate(file_id, "Output frequency", H5T_NATIVE_INT, space_one_id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(out_freq_dset_id >= 0);
    assert(H5Dwrite(out_freq_dset_id, H5T_NATIVE_INT, space_one_id, space_one_id, H5P_DEFAULT, &OUT_EVERY) >= 0);
    assert(H5Dclose(out_freq_dset_id) >= 0);
    assert(H5Sclose(space_one_id) >= 0);


    // ***** Coordinates *****
    constexpr hsize_t xdims = NX;
    constexpr hsize_t ydims = NY;

    const auto space_x_id = H5Screate_simple(1, &xdims, nullptr);
    const auto space_y_id = H5Screate_simple(1, &ydims, nullptr);
    assert(space_x_id >= 0);
    assert(space_y_id >= 0);

    const auto coords_group_id = H5Gcreate(file_id, "Coordinates",
                                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(coords_group_id >= 0); 

    const auto x_dset_id = H5Dcreate(coords_group_id, "x", H5T_NATIVE_DOUBLE, space_x_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto y_dset_id = H5Dcreate(coords_group_id, "y", H5T_NATIVE_DOUBLE, space_y_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(x_dset_id >= 0);
    assert(y_dset_id >= 0);

    assert(H5Dwrite(x_dset_id, H5T_NATIVE_DOUBLE, space_x_id, space_x_id, H5P_DEFAULT, x.data()) >= 0);
    assert(H5Dwrite(y_dset_id, H5T_NATIVE_DOUBLE, space_y_id, space_y_id, H5P_DEFAULT, y.data()) >= 0);

    assert(H5Dclose(x_dset_id)       >= 0);
    assert(H5Dclose(y_dset_id)       >= 0);
    assert(H5Gclose(coords_group_id) >= 0);


    // ***** Solution *****
    const auto u_group_id = H5Gcreate(file_id, "Solution",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u_group_id >= 0); 

    constexpr array<hsize_t, 2> udims{NX, NY};
    const auto space_u_id = H5Screate_simple(2, udims.data(), nullptr);
    assert(space_u_id >= 0);

    const auto u0_dset_id = H5Dcreate(u_group_id, "Iteration 0",
                                      H5T_NATIVE_DOUBLE, space_u_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u0_dset_id >= 0);
    assert(H5Dwrite(u0_dset_id, H5T_NATIVE_DOUBLE, space_u_id, space_u_id, H5P_DEFAULT, u.data()) >= 0);
    assert(H5Dclose(u0_dset_id) >= 0);



    /* ====================
     * Begin time evolution
     * ==================== */
    constexpr auto nxm1    = NX - 1;
    constexpr auto nym1    = NY - 1;
    constexpr auto nxm1_ny = ntot - NY;  // (NX - 1)*NY;

    constexpr double rx = ALPHA*DT/(dx*dx);
    constexpr double ry = ALPHA*DT/(dy*dy);

    constexpr double one_p_2rxpry = 1.0 + 2.0*(rx + ry);

    constexpr double rx_ux0 = rx*UX0;
    constexpr double rx_uxL = rx*UXL;
    constexpr double ry_uy0 = ry*UY0;
    constexpr double ry_uyL = ry*UYL;

    vector<double> u_old(ntot), r(ntot), p(ntot), Ap(ntot);

    constexpr hsize_t tdims = NT/OUT_EVERY;  // NOTE: integer division
    vector<double> time(tdims);

    // n == 0 has been written to file already
    for (auto n = decltype(NT){1}; n < NT; ++n) {
        u_old = u;

        for (auto i = decltype(NX){0}; i < NX; ++i) {
            const auto idx_i = i*NY;
            u_old[idx_i /*+0*/] += ry_uy0;  // u_old[i, 0]    += ry*UY0
            u_old[idx_i + nym1] += ry_uyL;  // u_old[i, NY-1] += ry*UYL
        }

        for (auto j = decltype(NY){0}; j < NY; ++j) {
            u_old[/* 0 + */ j] += rx_ux0;  // u_old[0,    j] += rx*UX0
            u_old[nxm1_ny + j] += rx_uxL;  // u_old[NX-1, j] += rx*UXL
        }


        /* Build the initial residual r = u_old - Au, with
         *     / 1+2(rx+ry)   -ry        0          0   ... \
         *     |    -rx    1+2(rx+ry)   -ry
         * A = |      0       -rx    1+2(rx+ry)   -ry   ... |
         *     |      .         .         .         .   ... |
         *     |      .         .         .         .   ... |
         *                                                                      */
        for (auto i = decltype(NX){0}; i < NX; ++i) {
            const auto idx_i  = i*NY;
            const auto idx_im = i*NY - NY;  // (i-1)*NY
            const auto idx_ip = i*NY + NY;  // (i+1)*NY

            for (auto j = decltype(NY){0}; j < NY; ++j) {
                const auto ij  = idx_i  + j;
                const auto imj = idx_im + j;
                const auto ipj = idx_ip + j;
                const auto ijm = idx_i  + j-1;
                const auto ijp = idx_i  + j+1;

                /* Build
                 *   r[ij] = u_old[ij] + rx*(u[imj] + u[ipj]) + ry*(u[ijm] + u[ijp]) - one_p_2rxpry*u[ij]
                 * at interior (i,j) points or suitable modifications of it
                 * otherwise                                                    */
                r[ij] = u_old[ij] - one_p_2rxpry*u[ij];

                if (i > 0) {
                    r[ij] += rx*u[imj];
                }

                if (i < nxm1) {
                    r[ij] += rx*u[ipj];
                }

                if (j > 0) {
                    r[ij] += ry*u[ijm];
                }

                if (j < nym1) {
                    r[ij] += ry*u[ijp];
                }
            }
        }

        p = r;


        /* There are at most NX*NY orthogonal directions in an NX*NY-dimensional
         * vector space                                                         */
        auto iter = decltype(ntot){0};

        for (iter = 0; iter < ntot; ++iter) {
            double rr_old = 0.0;
            double uu_old = 0.0;

            for (auto ij = decltype(ntot){0}; ij < ntot; ++ij) {
                rr_old +=     r[ij]*r[ij];
                uu_old += u_old[ij]*u_old[ij];
            }

            // Problem-scale-independent stopping criterion
            const auto r_rel = sqrt(rr_old/uu_old);

            if constexpr (VERYVERBOSE) {
                cout << "Conjugate gradient iteration " << iter
                     << "\tRelative residual: " << r_rel
                     << "\tTarget: " << CGTOL << endl;
            }

            if (r_rel < CGTOL) {
                break;
            }


            // Build A*p
            fill(Ap.begin(), Ap.end(), 0.0);

            for (auto i = decltype(NX){0}; i < NX; ++i) {
                const auto idx_i  = i*NY;
                const auto idx_im = i*NY - NY;  // (i-1)*NY
                const auto idx_ip = i*NY + NY;  // (i+1)*NY

                for (auto j = decltype(NY){0}; j < NY; ++j) {
                    const auto ij  = idx_i  + j;
                    const auto imj = idx_im + j;
                    const auto ipj = idx_ip + j;
                    const auto ijm = idx_i  + j-1;
                    const auto ijp = idx_i  + j+1;

                    Ap[ij] = one_p_2rxpry*p[ij];

                    if (i > 0) {
                        Ap[ij] -= rx*p[imj];
                    }

                    if (i < nxm1) {
                        Ap[ij] -= rx*p[ipj];
                    }

                    if (j > 0) {
                        Ap[ij] -= ry*p[ijm];
                    }

                    if (j < nym1) {
                        Ap[ij] -= ry*p[ijp];
                    }
                }
            }


            double pAp = 0.0;

            for (auto ij = decltype(ntot){0}; ij < ntot; ++ij) {
                pAp += p[ij]*Ap[ij];
            }

            const auto alpha = rr_old/pAp;

            for (auto ij = decltype(ntot){0}; ij < ntot; ++ij) {
                u[ij] += alpha*p[ij];
                r[ij] -= alpha*Ap[ij];
            }

            double rr = 0.0;

            for (auto ij = decltype(ntot){0}; ij < ntot; ++ij) {
                rr += r[ij]*r[ij];
            }

            const auto beta = rr/rr_old;

            for (auto ij = decltype(ntot){0}; ij < ntot; ++ij) {
                p[ij] = r[ij] + beta*p[ij];
            }
        }

        // If we got to this point without a solution u, we have a problem
        if (iter >= ntot) {  // > should never occur, but extra safety doesn't hurt
            ostringstream err_msg;
            err_msg << "Time " << n*DT
                    << ": couldn't find a solution to the linear system with the conjugate gradient method";
            throw runtime_error(err_msg.str());
            return 1;  // Not reached
        }


        // Write data to file
        static_assert(sizeof(std::size_t) <= sizeof(long long),
                      "size_t too big for long long");
        const auto [quot, rem] = lldiv(static_cast<long long>(n), static_cast<long long>(OUT_EVERY));
        if (rem == 0) {
            time[quot] = n*DT;

            ostringstream iteration_ss;
            iteration_ss << "Iteration " << n;

            const auto un_dset_id = H5Dcreate(u_group_id, iteration_ss.str().c_str(),
                                              H5T_NATIVE_DOUBLE, space_u_id,
                                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert(un_dset_id >= 0);
            assert(H5Dwrite(un_dset_id, H5T_NATIVE_DOUBLE, space_u_id, space_u_id, H5P_DEFAULT, u.data()) >= 0);
            assert(H5Dclose(un_dset_id) >= 0);

            if constexpr (VERBOSE) {
                cout << "Iteration " << n << ": data written to file" << endl;
            }
        }

        if constexpr (VERBOSE) {
            cout << "Iteration " << n << ": done" << endl;
        }
    }


    // Write time to file 
    const auto tspace_id = H5Screate_simple(1, &tdims, nullptr);
    assert(tspace_id >= 0);

    const auto t_dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, tspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(t_dset_id >= 0);

    assert(H5Dwrite(t_dset_id, H5T_NATIVE_DOUBLE, tspace_id, tspace_id, H5P_DEFAULT, time.data()) >= 0);

    assert(H5Dclose(t_dset_id) >= 0);
    assert(H5Sclose(tspace_id) >= 0);
 

    // Close up
    assert(H5Sclose(space_u_id) >= 0);
    assert(H5Gclose(u_group_id) >= 0);
    assert(H5Fclose(file_id) >= 0);

    cout << "Data written to '" << FILENAME << "'" << endl;

    return 0;
}
