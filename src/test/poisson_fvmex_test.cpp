#include "cfd/debug/printer.hpp"
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

using namespace cfd;

namespace {

class IProblem {
public:
    virtual ~IProblem() = default;
    virtual double exact_solution(Point p) const = 0;
    virtual double f(Point p) const = 0;

    // returns norm2(u - u_exact)
    double solve() {
        // 1. build SLAE
        CsrMatrix mat = approximate_lhs();
        std::vector<double> rhs = approximate_rhs();

        // 2. solve SLAE
        AmgcMatrixSolver solver; // iterative SLAE solver
        solver.set_matrix(mat);
        solver.solve(rhs, u_);

        // 3. compute norm2
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
            exact[i] = exact_solution(grid_->cell_center(i));
        }
        VtkUtils::add_cell_data(exact, "exact", filename);
    }

    std::shared_ptr<const IGrid> grid() const {
        return grid_;
    }

protected:
    IProblem(std::shared_ptr<IGrid> grid) : grid_(grid), ecol_(*grid_), u_(ecol_.size(), 0.0){};

    std::shared_ptr<IGrid> grid_;
    const FvmExtendedCollocations ecol_;
    std::vector<double> u_;

    virtual CsrMatrix approximate_lhs() const {
        LodMatrix mat(ecol_.size());

        // internal
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            // neighbour collocation points for the given face
            auto [icolloc, jcolloc] = ecol_.tab_face_colloc(iface);
            Point xi = ecol_.points[icolloc];
            Point xj = ecol_.points[jcolloc];
            double h = vector_abs(xj - xi);
            double v = grid_->face_area(iface) / h;
            mat.add_value(icolloc, icolloc, v);
            mat.add_value(jcolloc, jcolloc, v);
            mat.add_value(icolloc, jcolloc, -v);
            mat.add_value(jcolloc, icolloc, -v);
        }
        // boundary: dirichlet
        for (size_t icolloc: ecol_.face_collocations) {
            mat.set_unit_row(icolloc);
        }
        return mat.to_csr();
    }

    virtual std::vector<double> approximate_rhs() const {
        std::vector<double> rhs(ecol_.size(), 0.0);
        // internal
        for (size_t icolloc: ecol_.cell_collocations) {
            double value = f(ecol_.points[icolloc]);
            size_t icell = ecol_.cell_index(icolloc);
            double volume = grid_->cell_volume(icell);
            rhs[icell] = value * volume;
        }
        // boundary: dirichlet
        for (size_t iface: grid_->boundary_faces()) {
            // linear index of boundary collocation point
            size_t icolloc = ecol_.boundary_colloc(iface);
            rhs[icolloc] = exact_solution(grid_->face_center(iface));
        }
        return rhs;
    }

    double compute_norm2() const {
        double norm2 = 0;
        double full_area = 0;
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double diff = u_[icell] - exact_solution(grid_->cell_center(icell));
            norm2 += grid_->cell_volume(icell) * diff * diff;
            full_area += grid_->cell_volume(icell);
        }
        return std::sqrt(norm2 / full_area);
    }
};

/////////////////////////////////////////////////
// 2D Worker
/////////////////////////////////////////////////
class Problem2D : public IProblem {
public:
    Problem2D(std::shared_ptr<IGrid2D> grid) : IProblem(grid) {}

    // u(x) = sin(10*x^2)
    double exact_solution(Point p) const override {
        double x = p.x;
        double y = p.y;
        return sin(2 * M_PI * x + 3) + cos(10 * y * y);
    }

    // -d^2 u(x)/d x^2
    double f(Point p) const override {
        double x = p.x;
        double y = p.y;
        return 20.0 * sin(10 * y * y) + 400 * y * y * cos(10 * y * y) + 4 * M_PI * M_PI * sin(2 * M_PI * x + 3);
    }
};
} // namespace

TEST_CASE("Poisson 2D solver, Finite Volume Method", "[poisson2-fvmex]") {
    std::cout << std::endl << "--- [poisson2-fvmex] --- " << std::endl;

    auto grid = UnstructuredGrid2D::vtk_read_p(test_directory_file("tetragrid_500.vtk"), true);
    Problem2D worker(grid);

    // solve and find norm2
    double n2 = worker.solve();

    // save into poisson1_ncells={n_cells}.vtk
    worker.save_vtk("poisson2_fvmex.vtk");

    // CHECK
    CHECK(n2 == Approx(0.0633017).margin(1e-6));
}
