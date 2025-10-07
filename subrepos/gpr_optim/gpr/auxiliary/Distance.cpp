/*
 * Distance.cpp
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "Distance.h"

#include <cfloat>
#include <cmath>
#include <map>

#ifdef USE_HIGHS
#include "Highs.h"
#endif

namespace aux {

std::pair<gpr::Coord, gpr::Field<gpr::Index_t>> get_canonical_configuration(
    const gpr::Coord& x_in, const gpr::Field<gpr::Index_t>& types_in)
{
    const gpr::Index_t num_atoms = types_in.getSize();
    if (num_atoms == 0) {
        return {gpr::Coord(), gpr::Field<gpr::Index_t>()};
    }

    std::vector<AtomSignature> atom_signatures(num_atoms);

    // For each atom, compute its signature
    for (gpr::Index_t i = 0; i < num_atoms; ++i) {
        atom_signatures[i].original_index = i;
        atom_signatures[i].type = types_in(0, i);

        // The signature is the vector of inverse distances to all other atoms
        for (gpr::Index_t j = 0; j < num_atoms; ++j) {
            if (i == j) continue;
            double r_ij = (x_in.at(0, i) - x_in.at(0, j)).length();
            // Avoid division by zero, though unlikely in real systems
            if (r_ij > 1e-9) {
                atom_signatures[i].signature_vec.push_back(1.0 / r_ij);
            }
        }
        // Sort the signature vector itself to make it invariant to the ordering
        // of other atoms
        std::sort(atom_signatures[i].signature_vec.begin(),
                  atom_signatures[i].signature_vec.end());
    }

    // Sort the atoms based on their type and signature
    std::sort(atom_signatures.begin(), atom_signatures.end());

    // Create the new canonical configuration based on the sorted order
    gpr::Coord x_canonical;
    x_canonical.resize(1, num_atoms * 3);
    gpr::Field<gpr::Index_t> types_canonical(1, num_atoms);

    for (gpr::Index_t i = 0; i < num_atoms; ++i) {
        gpr::Index_t old_index = atom_signatures[i].original_index;
        x_canonical.set(0, i, x_in.at(0, old_index));
        types_canonical(0, i) = types_in(0, old_index);
    }

    return {x_canonical, types_canonical};
}

double solve_assignment_problem(const Eigen::MatrixXd& cost_matrix)
{
    const int n = cost_matrix.rows();
    if (n == 0) return 0.0;

    // We want to find a maximum weight matching, so we use max costs.
    // EMD is a minimum cost problem, so we negate the distances later.
    const Eigen::MatrixXd costs = -cost_matrix;

    Eigen::VectorXd lx = costs.rowwise().maxCoeff();
    Eigen::VectorXd ly = Eigen::VectorXd::Zero(n);
    std::vector<int> xy(n, -1);
    std::vector<int> yx(n, -1);
    std::vector<bool> S(n), T(n);

    auto match = [&](int u, auto& self) -> bool {
        S[u] = true;
        for (int v = 0; v < n; ++v) {
            if (!T[v] && std::abs(lx(u) + ly(v) - costs(u, v)) < 1e-9) {
                T[v] = true;
                if (yx[v] == -1 || self(yx[v], self)) {
                    yx[v] = u;
                    xy[u] = v;
                    return true;
                }
            }
        }
        return false;
    };

    for (int i = 0; i < n; ++i) {
        while (true) {
            std::fill(S.begin(), S.end(), false);
            std::fill(T.begin(), T.end(), false);
            if (match(i, match)) break;

            double delta = std::numeric_limits<double>::max();
            for (int u = 0; u < n; ++u) {
                if (S[u]) {
                    for (int v = 0; v < n; ++v) {
                        if (!T[v]) {
                            delta =
                                std::min(delta, lx(u) + ly(v) - costs(u, v));
                        }
                    }
                }
            }

            for (int j = 0; j < n; ++j) {
                if (S[j]) lx(j) -= delta;
                if (T[j]) ly(j) += delta;
            }
        }
    }

    double min_cost = 0.0;
    for (int i = 0; i < n; ++i) {
        min_cost += cost_matrix(i, xy[i]);
    }
    return min_cost;
}

#ifdef USE_HIGHS
// NOTE(rg):: Weirdly this is slower and spikier than the basic implementation
// above
double solve_assignment_problem_highs(const Eigen::MatrixXd& cost_matrix)
{
    const int n = cost_matrix.rows();
    if (n == 0) return 0.0;

    // 1. Create a HiGHS instance
    Highs highs;
    highs.setOptionValue("output_flag", false);  // Suppress solver output

    // 2. Define the LP model
    HighsModel model;
    model.lp_.num_col_ = n * n;  // n*n variables x_ij
    model.lp_.num_row_ = 2 * n;  // 2*n constraints

    // 3. Define the objective function (costs)
    // The variables are flattened: x_00, x_01, ..., x_0(n-1), x_10, ...
    model.lp_.col_cost_.resize(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            model.lp_.col_cost_[i * n + j] = cost_matrix(i, j);
        }
    }

    // 4. Define variable bounds (0 <= x_ij <= 1)
    model.lp_.col_lower_.assign(n * n, 0.0);
    model.lp_.col_upper_.assign(n * n, 1.0);

    // 5. Define the constraints
    // All constraint right-hand-sides are 1.
    model.lp_.row_lower_.assign(2 * n, 1.0);
    model.lp_.row_upper_.assign(2 * n, 1.0);

    // The constraint matrix A is very sparse. We define it in column-major
    // format. Each variable x_ij appears in exactly two constraints: row i and
    // row n+j.
    model.lp_.a_matrix_.start_.resize(n * n + 1);
    model.lp_.a_matrix_.index_.resize(2 * n * n);
    model.lp_.a_matrix_.value_.resize(2 * n * n);

    int nz = 0;  // Non-zero counter
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int var_idx = i * n + j;
            model.lp_.a_matrix_.start_[var_idx] = nz;

            // Entry for the "sum over j" constraint (row i)
            model.lp_.a_matrix_.index_[nz] = i;
            model.lp_.a_matrix_.value_[nz] = 1.0;
            nz++;

            // Entry for the "sum over i" constraint (row n+j)
            model.lp_.a_matrix_.index_[nz] = n + j;
            model.lp_.a_matrix_.value_[nz] = 1.0;
            nz++;
        }
    }
    model.lp_.a_matrix_.start_[n * n] = nz;
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;  // Specify format

    // 6. Solve the model
    highs.passModel(model);
    HighsStatus status = highs.run();

    if (status != HighsStatus::kOk) {
        // Handle error
        return std::numeric_limits<double>::max();
    }

    // 7. Return the optimal objective value
    return highs.getInfo().objective_function_value;
}
#endif

double calculate_aligned_rmsd(const gpr::Coord& P, const gpr::Coord& Q)
{
    // Assumes P and Q are single configurations (e.g., 1xN*3 vectors or Nx3
    // matrices)
    const gpr::Index_t num_atoms = P.getNumCols() / 3;
    if (num_atoms != Q.getNumCols() / 3) {
        // Handle error: structures must have the same number of atoms
        return std::numeric_limits<double>::max();
    }

    // Convert to a more convenient format, e.g., Eigen::MatrixX3d
    Eigen::MatrixX3d P_mat(num_atoms, 3);
    Eigen::MatrixX3d Q_mat(num_atoms, 3);
    for (gpr::Index_t i = 0; i < num_atoms; ++i) {
        P_mat.row(i) << P(0, 3 * i), P(0, 3 * i + 1), P(0, 3 * i + 2);
        Q_mat.row(i) << Q(0, 3 * i), Q(0, 3 * i + 1), Q(0, 3 * i + 2);
    }

    // 1. Center the coordinates
    Eigen::RowVector3d p_centroid = P_mat.colwise().mean();
    Eigen::RowVector3d q_centroid = Q_mat.colwise().mean();
    P_mat.rowwise() -= p_centroid;
    Q_mat.rowwise() -= q_centroid;

    // 2. Calculate the covariance matrix H = P_centered^T * Q_centered
    Eigen::Matrix3d H = P_mat.transpose() * Q_mat;

    // 3. Compute SVD of H
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d V = svd.matrixV();
    Eigen::Matrix3d U = svd.matrixU();

    // 4. Calculate the optimal rotation matrix R
    double d = (V * U.transpose()).determinant();
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    D(2, 2) = d;
    Eigen::Matrix3d R = V * D * U.transpose();

    // 5. Apply rotation to P and calculate RMSD
    P_mat = (R * P_mat.transpose()).transpose();

    double sum_sq_dist = (P_mat - Q_mat).squaredNorm();
    return std::sqrt(sum_sq_dist / num_atoms);
}

Distance::Distance() { }

Distance::~Distance() { }

/**
 * @brief Calculates the aligned Root Mean Square Deviation (RMSD) between
 * sets of configurations.
 */
void Distance::dist_rmsd(const gpr::Coord& x1, const gpr::Coord& x2,
                         const gpr::AtomsConfiguration& conf_info,
                         gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();  // Number of configurations in x1
    gpr::Index_t n2 = x2.getNumRows();  // Number of configurations in x2
    gpr::Index_t num_cols = x1.getNumCols();

    // Basic check for consistent atom numbers
    if (num_cols != x2.getNumCols()) {
        // Or throw an exception for a more robust error handling
        dist.clear();
        return;
    }

    dist.resize(n1, n2);
    const gpr::Field<gpr::Index_t>& types = conf_info.atoms_mov.type;

    for (gpr::Index_t n = 0; n < n1; ++n) {
        // Manually extract the n-th configuration (row) from x1 into a new
        // single-row Coord object 'p'.
        gpr::Coord p;
        p.resize(1, num_cols);
        for (gpr::Index_t k = 0; k < num_cols; ++k) {
            p(0, k) = x1(n, k);
        }

        auto p_canonical_pair = get_canonical_configuration(p, types);

        for (gpr::Index_t m = 0; m < n2; ++m) {
            // Manually extract the m-th configuration (row) from x2 into 'q'.
            gpr::Coord q;
            q.resize(1, num_cols);
            for (gpr::Index_t k = 0; k < num_cols; ++k) {
                q(0, k) = x2(m, k);
            }

            auto q_canonical_pair = get_canonical_configuration(q, types);

            // Calculate the aligned RMSD and store it in the distance matrix.
            dist(n, m) = calculate_aligned_rmsd(p_canonical_pair.first,
                                                q_canonical_pair.first);
        }
    }
}

void Distance::dist_max1Dlog(const gpr::Coord& x1, const gpr::Coord& x2,
                             const gpr::AtomsConfiguration& conf_info,
                             gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = x1.getNumCols() / 3;
    gpr::Index_t N_fro = conf_info.atoms_froz_active.positions.getNumCols() / 3;

    dist.resize(n1, n2);

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t i = 0; i < N_mov - 1; ++i) {
            for (gpr::Index_t j = i + 1; j < N_mov; ++j) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij_1 = (x1.at(n, j) - x1.at(n, i)).length();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double r_ij_2 = (x2.at(m, j) - x2.at(m, i)).length();
                        double tmp = fabs(log(r_ij_2 / r_ij_1));
                        if (tmp > dist(n, m)) dist(n, m) = tmp;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .length();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double r_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .length();
                        double tmp = fabs(log(r_ij_2 / r_ij_1));
                        if (tmp > dist(n, m)) dist(n, m) = tmp;
                    }
                }
            }
        }
    }
}

void Distance::dist_at(const gpr::Coord& x1, const gpr::Coord& x2,
                       const gpr::AtomsConfiguration& conf_info,
                       const gpr::Field<double>& lengthscale,
                       gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getSize();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getSize();
    gpr::Field<double> s2;
    gpr::Index_t size_s2 = lengthscale.getSize();

    dist.resize(n1, n2);

    s2.resize(1, size_s2);
    for (gpr::Index_t n = 0; n < s2.getSize(); ++n)
        s2[n] = 1. / (lengthscale[n] * lengthscale[n]);

    // If ARD is not used make s a vector of equal elements
    if (size_s2 == 1) {
        double ref_s2 = s2(0, 0);
        size_s2 = conf_info.n_pt;
        s2.resize(1, size_s2);
        for (gpr::Index_t n = 0; n < size_s2; ++n) {
            s2(0, n) = ref_s2;
        }
    }

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t j = 0; j < N_mov - 1; ++j) {
            for (gpr::Index_t i = j + 1; i < N_mov; ++i) {
                double s2_val =
                    s2(0, conf_info.pairtype(conf_info.atoms_mov.type(0, i),
                                             conf_info.atoms_mov.type(0, j)));
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 = (x1.at(n, j) - x1.at(n, i)).rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) - x2.at(m, i)).rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist(n, m) += 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                double s2_val = s2(
                    0,
                    conf_info.pairtype(conf_info.atoms_froz_active.type(0, i),
                                       conf_info.atoms_mov.type(0, j)));
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist(n, m) += 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    for (gpr::Index_t n = 0; n < dist.getSize(); ++n) {
        dist[n] = sqrt(dist[n]);
    }
}

void Distance::dist_at_vec(const gpr::Coord& x1, const gpr::Coord& x2,
                           const gpr::AtomsConfiguration& conf_info,
                           const gpr::Field<double>& lengthscale,
                           std::vector<gpr::Field<double>>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getSize();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getSize();
    gpr::Field<double> s2;
    gpr::Index_t size_s2 = lengthscale.getSize();

    dist.resize(conf_info.n_pt);
    for (gpr::Index_t n = 0; n < dist.size(); ++n)
        dist[n].resize(n1, n2);

    s2.resize(1, size_s2);
    for (gpr::Index_t n = 0; n < size_s2; ++n)
        s2(0, n) = 1. / (lengthscale(0, n) * lengthscale(0, n));

    // If ARD is not used make s a vector of equal elements
    if (size_s2 == 1) {
        double ref_s2 = s2(0, 0);
        size_s2 = conf_info.n_pt;
        s2.resize(1, size_s2);
        for (gpr::Index_t n = 0; n < size_s2; ++n) {
            s2(0, n) = ref_s2;
        }
    }

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t j = 0; j < N_mov - 1; ++j) {
            for (gpr::Index_t i = j + 1; i < N_mov; ++i) {
                gpr::Index_t pt =
                    conf_info.pairtype(conf_info.atoms_mov.type(0, i),
                                       conf_info.atoms_mov.type(0, j));
                double s2_val = s2(0, pt);
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 = (x1.at(n, j) - x1.at(n, i)).rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) - x2.at(m, i)).rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist[pt](n, m) -= 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                gpr::Index_t pt =
                    conf_info.pairtype(conf_info.atoms_froz_active.type(0, i),
                                       conf_info.atoms_mov.type(0, j));
                double s2_val = s2(0, pt);
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist[pt](n, m) -= 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    //    for(gpr::Index_t n = 0; n < n1; ++n) {
    //        for(gpr::Index_t m = 0; m < n2; ++m) {
    //            dist(n, m) = sqrt(dist(n, m));
    //        }
    //    }
}

void Distance::mindist_interatomic(const gpr::Coord& x,
                                   const gpr::AtomsConfiguration& conf_info,
                                   gpr::Field<double>& dist)
{
    // Use clearer variable names and correct size calculations
    const auto num_configs = x.getNumRows();
    const auto num_moving_atoms = x.getNumCols() / 3;
    const auto& frozen_pos = conf_info.atoms_froz_active.positions;
    const auto num_frozen_atoms = frozen_pos.getNumCols() / 3;

    // Handle the edge case of no moving atoms
    if (num_moving_atoms == 0) {
        dist.clear();
        return;
    }

    dist.resize(num_configs, num_moving_atoms);
    dist.set(DBL_MAX);

    // Loop over configurations first to improve memory access patterns
    // (cache-friendlier)
    for (gpr::Index_t n = 0; n < num_configs; ++n) {
        // --- Distances between moving atoms for configuration 'n' ---
        if (num_moving_atoms > 1) {
            for (gpr::Index_t i = 0; i < num_moving_atoms - 1; ++i) {
                for (gpr::Index_t j = i + 1; j < num_moving_atoms; ++j) {
                    const double r_ij = (x.at(n, j) - x.at(n, i)).length();

                    // Update the minimum distance for both atoms in the pair
                    // Use std::min for clarity and correctness.
                    dist(n, i) = std::min(dist(n, i), r_ij);
                    dist(n, j) = std::min(dist(n, j), r_ij);
                }
            }
        }

        // --- Distances from moving to frozen atoms for configuration 'n' ---
        if (num_frozen_atoms > 0) {
            for (gpr::Index_t j = 0; j < num_moving_atoms;
                 ++j) {  // For each moving atom
                for (gpr::Index_t i = 0; i < num_frozen_atoms;
                     ++i) {  // Check against each frozen atom
                    const double r_ij =
                        (x.at(n, j) - frozen_pos.at(0, i)).length();

                    // Update the minimum distance ONLY for the current moving
                    // atom 'j'
                    dist(n, j) = std::min(dist(n, j), r_ij);
                }
            }
        }
    }
}

void Distance::dist_emd(const gpr::Coord& x1, const gpr::Coord& x2,
                        const gpr::AtomsConfiguration& conf_info,
                        gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();  // Number of configurations in x1
    gpr::Index_t n2 = x2.getNumRows();  // Number of configurations in x2
    gpr::Index_t num_atoms = x1.getNumCols() / 3;

    if (num_atoms != x2.getNumCols() / 3) {
        dist.clear();
        return;
    }

    dist.resize(n1, n2);
    const gpr::Field<gpr::Index_t>& types = conf_info.atoms_mov.type;

    // Group atom indices by type
    std::map<gpr::Index_t, std::vector<gpr::Index_t>> indices_by_type;
    for (gpr::Index_t i = 0; i < types.getSize(); ++i) {
        indices_by_type[types(0, i)].push_back(i);
    }

    for (gpr::Index_t n = 0; n < n1; ++n) {
        for (gpr::Index_t m = 0; m < n2; ++m) {
            double max_of_avg_type_emds = 0.0;
            double total_emd_for_type = 0.0;

            // Calculate EMD for each atom type separately
            for (auto const& [type, indices]: indices_by_type) {
                int num_type_atoms = indices.size();
                if (num_type_atoms == 0) continue;

                // Build the cost matrix for this type
                Eigen::MatrixXd cost_matrix(num_type_atoms, num_type_atoms);
                for (int i = 0; i < num_type_atoms; ++i) {
                    for (int j = 0; j < num_type_atoms; ++j) {
                        gpr::Index_t atom_idx1 = indices[i];
                        gpr::Index_t atom_idx2 = indices[j];
                        cost_matrix(i, j) =
                            (x1.at(n, atom_idx1) - x2.at(m, atom_idx2))
                                .length();
                    }
                }

                // Solve the assignment problem for this atom type
#ifdef USE_HIGHS
                total_emd_for_type =
                    solve_assignment_problem_highs(cost_matrix);
#else
                total_emd_for_type = solve_assignment_problem(cost_matrix);
#endif
                // Calculate the AVERAGE displacement for THIS type
                double avg_emd_for_type = total_emd_for_type / num_type_atoms;

                // Keep track of the MAXIMUM average we've seen
                if (avg_emd_for_type > max_of_avg_type_emds) {
                    max_of_avg_type_emds = avg_emd_for_type;
                }
            }

            // NOTE(rg): This is now INTENSIVE, so on average no atom group
            // should move more than this amount
            dist(n, m) = max_of_avg_type_emds;
        }
    }
}

} /* namespace aux */
