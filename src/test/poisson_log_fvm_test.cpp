#include "cfd/debug/printer.hpp"
#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd27_test.hpp"
#include "utils/filesystem.hpp"
#include <ranges>

using namespace cfd;

namespace {

struct ITestData {
    virtual ~ITestData() = default;
    virtual double exact_solution(Point p) const = 0;
    virtual double lambda(Point p) const = 0;
    virtual double f(Point p) const = 0;
};

class Problem {
public:
    Problem(std::shared_ptr<IGrid> grid, std::shared_ptr<ITestData> test_data)
        : grid_{grid},
          test_data_{test_data},
          lambda_cell_{std::views::iota(0uz, grid_->n_cells()) | std::views::transform([&](size_t icell) {
                           return test_data->lambda(grid_->cell_center(icell));
                       }) |
                       std::ranges::to<std::vector<double>>()},
          lambda_face_{Problem::build_lambda_face(lambda_cell_, grid_)},
          ecol_(*grid_),
          grad_computer_(std::make_shared<LeastSquaresFvmFaceGradient>(*grid, ecol_)),
          u_(ecol_.size(), 0.0){};

    virtual ~Problem() = default;

    // returns norm2(u - u_exact)
    double solve(size_t n_iters) {
        // 1. build SLAE
        CsrMatrix mat = approximate_lhs();

        // 2. solve SLAE
        AmgcMatrixSolver solver; // iterative SLAE solver
        solver.set_matrix(mat);

        // 3. correction loop
        for (size_t iter = 0; iter < n_iters; ++iter) {
            std::vector<double> rhs = approximate_rhs_skew_corrected(u_);
            solver.solve(rhs, u_);
        }

        // 4. compute norm2
        return compute_norm2();
    }

    // saves numerical and exact solution into the vtk format
    void save_vtk(const std::string& filename) const {
        // save grid
        grid_->save_vtk(filename);

        // save numerical solution
        VtkUtils::add_cell_data(u_, "numerical", filename, grid_->n_cells());

        // save exact solution
        std::vector<double> exact(grid_->n_cells());
        for (size_t i = 0; i < grid_->n_cells(); ++i) {
            exact[i] = test_data_->exact_solution(grid_->cell_center(i));
        }
        VtkUtils::add_cell_data(exact, "exact", filename);
    }

    std::shared_ptr<const IGrid> grid() const {
        return grid_;
    }

protected:
    const std::shared_ptr<IGrid> grid_;
    const std::shared_ptr<ITestData> test_data_;
    const std::vector<double> lambda_cell_;
    const std::vector<double> lambda_face_;
    const FvmExtendedCollocations ecol_;
    const std::shared_ptr<IFvmFaceGradient> grad_computer_;
    std::vector<double> u_;

    static std::vector<double> build_lambda_face(const std::vector<double>& lambda_cell, std::shared_ptr<IGrid> grid) {
        // ------- compute lambda_face
        std::vector<double> lambda_face(grid->n_faces());
        // internal
        for (size_t iface: grid->internal_faces()) {
            auto [icell, jcell] = grid->tab_face_cell(iface);
            // arithmetic average value on face
            lambda_face[iface] = 0.5 * (lambda_cell[icell] + lambda_cell[jcell]);
        }
        // boundary
        for (size_t iface: grid->boundary_faces()) {
            auto [minus, plus] = grid->tab_face_cell(iface);
            size_t icell = (minus == INVALID_INDEX) ? plus : minus;
            lambda_face[iface] = lambda_cell[icell];
        }
        return lambda_face;
    }

    virtual CsrMatrix approximate_lhs() const {
        LodMatrix mat(ecol_.size());

        // all faces
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [icolloc, jcolloc] = ecol_.tab_face_colloc(iface);
            Point xi = ecol_.points[icolloc];
            Point xj = ecol_.points[jcolloc];
            double h = vector_abs(xj - xi);

            // skew method: rotate
            double cbs = 1.0;

            double v = cbs * lambda_face_[iface] * grid_->face_area(iface) / h;
            mat.add_value(icolloc, icolloc, v);
            mat.add_value(jcolloc, jcolloc, v);
            mat.add_value(icolloc, jcolloc, -v);
            mat.add_value(jcolloc, icolloc, -v);
        }

        // boundary fix: dirichlet
        for (size_t icolloc: ecol_.face_collocations) {
            mat.set_unit_row(icolloc);
        }
        return mat.to_csr();
    }

    virtual std::vector<double> approximate_rhs_skew_corrected(const std::vector<double>& u) const {
        // init
        std::vector<double> rhs(ecol_.size(), 0.0);

        // skew
        std::vector<Vector> gradu = grad_computer_->compute(u);
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [i, j] = ecol_.tab_face_colloc(iface);
            double area = grid_->face_area(iface);
            Vector cij = ecol_.points[j] - ecol_.points[i];

            // skew method: rotate
            Vector c = cij / vector_abs(cij);

            Vector a = grid_->face_normal(iface) - c;
            double val = dot_product(gradu[iface], a) * area;
            rhs[i] += val;
            rhs[j] -= val;
        }

        // boundary: dirichlet
        for (size_t iface: grid_->boundary_faces()) {
            // linear index of boundary collocation point
            size_t icolloc = ecol_.boundary_colloc(iface);
            rhs[icolloc] = test_data_->exact_solution(grid_->face_center(iface));
        }
        return rhs;
    }

    double compute_norm2() const {
        double norm2 = 0;
        double full_area = 0;
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double diff = u_[icell] - test_data_->exact_solution(grid_->cell_center(icell));
            norm2 += grid_->cell_volume(icell) * diff * diff;
            full_area += grid_->cell_volume(icell);
        }
        return std::sqrt(norm2 / full_area);
    }
};

// TEST data: logarithmic function in r=[0.05, 1.05]
struct TestDataRadial : public ITestData {
    double lambda(Point p) const override {
        double r = std::hypot(p.x, p.y);
        return (r < RMID) ? LAMBDA0 : LAMBDA1;
    }

    double exact_solution(Point p) const override {
        static constexpr double B = std::log(RMID / R0) / LAMBDA0 + std::log(R1 / RMID) / LAMBDA1;
        double r = std::hypot(p.x, p.y);
        double a = 1.0;
        double b = 0.0;
        if (r < RMID) {
            return a + (b - a) / (B * LAMBDA0) * std::log(r / R0);
        } else {
            return b + (b - a) / (B * LAMBDA1) * std::log(r / R1);
        }
    }

    double f(Point) const override {
        return 0.0;
    }

private:
    static constexpr double R0 = 0.05;
    static constexpr double R1 = 1.05;
    static constexpr double RMID = (R0 + R1) / 2.0;
    static constexpr double LAMBDA0 = 10.0;
    static constexpr double LAMBDA1 = 1.0;
};

} // namespace

TEST_CASE("Poisson 2D solver for radial problem, Finite Volume Method", "[poisson-log-fvm-2d]") {
    std::cout << std::endl << "--- [poisson-log-fvm-2d] --- " << std::endl;
    auto test_data = std::make_shared<TestDataRadial>();

    auto grid = UnstructuredGrid2D::vtk_read_p(test_directory_file("radial_500.vtk"), true);
    Problem worker(grid, test_data);

    // solve and find norm2
    double n2 = worker.solve(10);

    // save into vtk
    worker.save_vtk("poisson_log_fvm.vtk");

    // CHECK
    CHECK(n2 == Approx(0.0205646191).margin(1e-6));
}
